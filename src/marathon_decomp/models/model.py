"""The production marathon-decomposition model and its ALS fitter.

Model (centered form):

    log t_ij = u_i + v_j
             + aging(A^n_ij)                          # basis expansion
             + entry_age(A^e_i, A^n_ij)               # gamma block
             + d_i*(A^n_ij - mean_i(A^n))             # per-athlete drift
             + eps_ij,        eps_ij ~ t_nu(0, sigma^2)

    u_i  athlete ability        v_j  race difficulty (gauge: sum_j v_j = 0)
    A^n  career age (years since the athlete's first race in the slice)
    A^e  entry age (years between birth and first race)

The aging block is a basis expansion of A^n (see marathon_decomp.aging):

    aging(A_n) = sum_k theta_aging_k * B_k(A_n)

    basis_kind="poly":   {A_n, ..., A_n^degree}, optionally QR-orthogonalized
    basis_kind="spline": natural cubic regression spline

The entry-age (gamma) block has two forms:

    gamma_form="scalar":   gamma * (A^e_i - mean) * A^n_ij            (1 col)
    gamma_form="varying":  sum_k gamma_k * (A^e_i - mean) * B_k(A^n)  (K cols)

EB prior d_i ~ N(0, omega_d^2) over eligible athletes only, an optional
runner-reliability weight w_ath_i = n_i/(n_i+n0), and nu = inf recovers L2.

Algorithm: block coordinate descent with a running residual, IRLS for the
Student-t loss, a joint solve for the {aging, gamma} block, and an EB update
for omega_d^2. ModelAnderson (anderson.py) wraps the same per-iteration step
body with Anderson acceleration. The full derivation lives in
docs/model_description_v2.md.

Athletes with no DOB info (A_e = NaN) contribute zero to the entry-age block
at all their cells but participate in every other block.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
from scipy.special import gammaln

from ..aging import build_poly_basis, default_knots_from_An, ncs_basis
from ..base import BaseModel, FitResult
from ..data import FitData
from ..kernels import cd_inner_scalars, irls_weights

_EPS = 1e-12
_SEC_PER_YEAR = 365.25 * 86400.0   # same convention as data.py

# Default iter-0 value for the EB drift prior variance when neither omega_d2_init
# nor omega_d2_fixed is given. The converged omega_d2 is init-INVARIANT (the EM
# update is a contraction to the unimodal type-II MLE), so this affects only
# convergence speed/robustness, not results. Chosen ~just below the observed
# omega* ~ 3e-4 (stable across all production slices) and ~4 orders ABOVE the
# near-zero EM stall manifold: a too-small start (e.g. <~1e-7) crawls and can hit
# max_iter without converging. See scripts/02_model_selection/athlete_drift
# (e02_omega_init) for the init-sensitivity study behind this value.
_OMEGA_D2_INIT_DEFAULT = 1e-4


# ---------------------------------------------------------------------------
# Config dataclasses
# ---------------------------------------------------------------------------

@dataclass
class ModelConfig:
    """What the model *is*. Pure modeling knobs, no numerics."""

    # term toggles. use_phi12 toggles the whole aging-basis block; use_gamma
    # toggles the entry-age block; use_d toggles per-athlete drift.
    use_phi12: bool = True
    use_gamma: bool = True
    use_d: bool = True

    # --- aging basis (the phi_k block) ---
    basis_kind: Literal["poly", "spline"] = "poly"
    degree: int = 2               # poly only; degree=2 == classic phi1,phi2
    orthogonalize: bool = True    # poly only; QR-orthonormalize columns
    n_knots: int = 5              # spline only; used when `knots` is None
    knots: tuple[float, ...] | None = None   # spline only; explicit knots
    gamma_form: Literal["scalar", "varying"] = "scalar"
    # Aging-curve fit domain. None => the aging/gamma coefficients are fit on
    # all cells. A float caps the fit to cells with A_n <= aging_an_max: those
    # coefficients see only A_n<=cap, but the curve still predicts (and the
    # running residual still patches) every cell, so A_n>cap cells stay in the
    # rest of the estimation and the curve simply extrapolates beyond the cap.
    aging_an_max: float | None = None

    # eligibility for d_i
    d_min_n: int = 3                 # >=3 finishes needed for a slope
    d_min_span_years: float = 24/52  # >= two 12-week blocks of career span

    # EB prior on d_i. None => estimate omega_d2 via type-II MLE; a fixed value
    # freezes it. omega_d2_max caps the estimate (omega_d2 is the variance of a
    # per-year log-time drift: d_i ~ 0.05 means ~5% improvement per year).
    # omega_d2_init = iter-0 guess only (result is init-invariant; affects
    # convergence speed). None -> _OMEGA_D2_INIT_DEFAULT (=1e-4).
    omega_d2_init: float | None = None
    omega_d2_fixed: float | None = None
    omega_d2_max: float | None = None
    freeze_eb_prior: bool = False

    # noise model
    nu: float = float("inf")              # inf = Gaussian; finite = Student-t
    sigma2_init: float | None = None

    # runner-reliability weight. 0 = off.
    runner_weight_n0: float = 0.0

    # gauge
    gauge: Literal["mean_v_zero"] = "mean_v_zero"

    # APC gauge (§4.2). When True, enforce the G1 (β=0) gauge every iteration:
    # force the secular trend of v_j on race date to zero by shifting the linear
    # aging slope. Predictions are invariant; the gauge is re-applied post-fit
    # via `Model.apply_apc_gauge_beta0()` regardless of this flag.
    apc_gauge_beta0: bool = True


@dataclass
class FitterConfig:
    """How the model is *solved*. Pure numerics, no modeling content."""

    max_outer_iter: int = 300
    tol: float = 1e-10
    # Stopping scalar (both use a relative tolerance against `tol`):
    #   "loglik": |dl| < tol*(|l|+eps) — tracks the actual objective; correct
    #             under Student-t and EB updates that don't change RSS. Default.
    #   "rss":    |dRSS| < tol*(|RSS|+eps) — only valid in the Gaussian limit.
    stop_criterion: Literal["loglik", "rss"] = "loglik"
    cd_inner_iters: int = 8
    # Aging-block solver:
    #   "direct":   one dense np.linalg.solve on the (Kblk x Kblk) weighted
    #               Gram per outer iter — exact joint optimum. Default and
    #               required for spline bases (correlated columns trap CD).
    #   "numba_cd": cd_inner_scalars sweep; only valid for poly bases.
    inner_solver: Literal["numba_cd", "direct"] = "direct"
    residual_refresh_every: int = 1
    # "mean": u_i = mean response of athlete i, v_j = 0 (sensible sigma2 scale).
    # "zeros" / "random": for debugging / determinism.
    # "warmstart": load u, v, d, theta_aging, gamma from `warmstart`.
    init: Literal["mean", "zeros", "random", "warmstart"] = "mean"
    seed: int = 0
    warmstart: dict[str, np.ndarray] | None = None
    verbose: int = 0
    record_trace: bool = False


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------

class Model(BaseModel):
    """Baseline + aging + drift model, fit by block coordinate descent (ALS)."""

    NAME = "model_als"
    # 1.3.0: fixed the _t_j units bug that silently disabled the beta=0 APC
    # gauge (see _setup_aging / apply_apc_gauge_beta0). Fits <= 1.2.0 are
    # off-manifold; re-gauge their v post-hoc or refit with >= 1.3.0.
    VERSION = "1.3.0"
    CONFIG_CLS = ModelConfig
    PARAM_SPEC = {
        "u": ("I",), "v": ("J",), "d": ("I",),
        "theta_aging": ("K_basis",),
        "gamma": ("K_gamma",),     # length 1 for scalar, K_basis for varying
        "sigma2": (), "omega_d2": (),
        "nu": (),
    }

    # Optional per-cell resampling weight (length N), folded into the IRLS
    # weight so it flows into every weighted block. Used by the athlete
    # bootstrap (marathon_decomp.resample). None => no-op.
    boot_w_cell: np.ndarray | None = None

    def __init__(
        self,
        data: FitData,
        config: ModelConfig | None = None,
        fitter: FitterConfig | None = None,
    ):
        super().__init__(data, config or ModelConfig())
        self.fitter = fitter or FitterConfig()
        self._precompute()

    # ------------------------------------------------------------------
    # Precomputation (one-shot, data-only)
    # ------------------------------------------------------------------
    def _precompute(self) -> None:
        d = self.data
        cfg = self.config
        row, col = d.row_idx, d.col_idx
        I, J, N = d.I, d.J, d.N

        self.n_per_ath = np.bincount(row, minlength=I).astype(np.float64)
        self.n_per_race = np.bincount(col, minlength=J).astype(np.float64)

        # A^n centering for d_i: A_tilde^n_ij = A^n_ij - mean_i(A^n).
        sum_An = np.bincount(row, weights=d.A_n, minlength=I)
        An_bar = np.where(self.n_per_ath > 0, sum_An / np.maximum(self.n_per_ath, 1.0), 0.0)
        self.An_tilde_cell = d.A_n - An_bar[row]

        # A^e centering for gamma. Athletes with NaN A_e are excluded from the
        # mean and get (A^e - mean) = 0 (zero contribution to gamma, full
        # participation elsewhere).
        known = ~np.isnan(d.A_e)
        Ae_bar = float(d.A_e[known].mean()) if known.any() else 0.0
        Ae_c = np.where(known, d.A_e - Ae_bar, 0.0)
        self.Ae_bar = Ae_bar
        self.Ae_c_cell = Ae_c[row]
        self.has_ae_cell = known[row]

        # Eligibility for d_i: n_i >= d_min_n AND career span >= d_min_span.
        if cfg.use_d:
            An_min = np.full(I, np.inf)
            An_max = np.full(I, -np.inf)
            np.minimum.at(An_min, row, d.A_n)
            np.maximum.at(An_max, row, d.A_n)
            span = np.where(np.isfinite(An_max) & np.isfinite(An_min), An_max - An_min, 0.0)
            self.eligible_d = (
                (self.n_per_ath >= cfg.d_min_n)
                & (span >= cfg.d_min_span_years)
            )
        else:
            self.eligible_d = np.zeros(I, dtype=bool)

        # ---- aging basis B (N, K_basis) ------------------------------
        if cfg.basis_kind == "poly":
            B, transform = build_poly_basis(d.A_n, cfg.degree, cfg.orthogonalize)
            self._basis_transform = transform
            self._spline_knots = None
        elif cfg.basis_kind == "spline":
            if cfg.knots is not None:
                knots = np.asarray(cfg.knots, dtype=np.float64)
            else:
                knots = default_knots_from_An(d.A_n, n_knots=cfg.n_knots)
            B = ncs_basis(d.A_n, knots)
            self._basis_transform = None
            self._spline_knots = knots
        else:
            raise ValueError(f"unknown basis_kind={cfg.basis_kind!r}")

        if B.shape[0] != N:
            raise ValueError(f"basis row count {B.shape[0]} != N={N}")
        self.K_basis = int(B.shape[1])
        self._B_aging = B   # kept for curve reconstruction / varying-gamma cols

        # Aging-curve fit-domain mask (zeros the aging-block weight beyond cap).
        if cfg.aging_an_max is not None:
            self._aging_cell_w = (d.A_n <= float(cfg.aging_an_max)).astype(np.float64)
        else:
            self._aging_cell_w = None

        # ---- assemble x_scalar (K_blk, N) for the joint aging block --
        # rows [0:K_basis)     = B columns (the phi_k regressors)
        # rows [K_basis:K_blk) = gamma columns (1 or K_basis)
        if cfg.gamma_form == "scalar":
            K_gamma = 1
            gamma_cols = (self.Ae_c_cell * d.A_n)[None, :]            # (1, N)
        elif cfg.gamma_form == "varying":
            K_gamma = self.K_basis
            gamma_cols = self.Ae_c_cell[None, :] * B.T                # (K, N)
        else:
            raise ValueError(f"unknown gamma_form={cfg.gamma_form!r}")
        self.K_gamma = K_gamma
        self.K_blk = self.K_basis + K_gamma

        x_scalar = np.empty((self.K_blk, N), dtype=np.float64)
        x_scalar[: self.K_basis] = B.T
        x_scalar[self.K_basis :] = gamma_cols
        self.x_scalar = x_scalar

        mask = np.empty(self.K_blk, dtype=np.bool_)
        mask[: self.K_basis] = bool(cfg.use_phi12)
        mask[self.K_basis :] = bool(cfg.use_gamma)
        self.scalar_use_mask = mask
        # The mask is static, so cache the active regressor rows once
        # (contiguous) — the direct aging solve reuses this every iteration
        # instead of re-copying x_scalar[mask].
        self._x_active = np.ascontiguousarray(self.x_scalar[mask])

        # Runner-reliability weight.
        n0 = float(cfg.runner_weight_n0)
        w_ath = self.n_per_ath / (self.n_per_ath + n0) if n0 > 0 else np.ones(I)
        self.w_ath_cell = w_ath[row]

        # APC gauge quantities: race date and debut date in fractional Julian years,
        # satisfying t_j[col_idx] - f_i[row_idx] = A_n exactly.
        # BUG FIX (fitter 1.3.0): the previous line read
        #     t_j = d.race_date.astype("int64") * (1e-9 / _SEC_PER_YEAR)
        # which assumed race_date was datetime64[ns]. It is datetime64[s] (see
        # data.py), so .astype("int64") already yields SECONDS, and the extra
        # 1e-9 shrank t_j by ~1e9 (e.g. year 2025 -> ~5e-8 instead of ~55.6).
        # var(t_j) then fell below _EPS, so apply_apc_gauge_beta0()'s
        # degenerate-date guard returned early on every call -> the beta=0 APC
        # gauge SILENTLY NEVER RAN. All fits saved by fitter <= 1.2.0 therefore
        # sit OFF the beta=0 manifold (slope(v ~ t) != 0); their v carries a
        # prior-pinned date tilt that downstream scripts must remove post-hoc.
        # Cast to [s] first so the conversion is datetime-resolution-robust.
        t_j = d.race_date.astype("datetime64[s]").astype("int64") / _SEC_PER_YEAR
        f_i_cell = t_j[col] - d.A_n    # debut date in years per cell: t_race - A_n
        f_i_sum = np.bincount(row, weights=f_i_cell, minlength=I)
        f_i = np.where(self.n_per_ath > 0, f_i_sum / np.maximum(self.n_per_ath, 1.0), 0.0)
        self._t_j = t_j   # (J,) race dates in fractional Julian years
        self._f_i = f_i   # (I,) debut dates in fractional Julian years

    # ------------------------------------------------------------------
    # Initialization
    # ------------------------------------------------------------------
    def _init_params(self) -> None:
        d, cfg, fcfg = self.data, self.config, self.fitter
        I, J, N = d.I, d.J, d.N
        rng = np.random.default_rng(fcfg.seed)

        if fcfg.init == "warmstart" and fcfg.warmstart is not None:
            ws = fcfg.warmstart
            u = ws.get("u", np.zeros(I)).astype(np.float64).copy()
            v = ws.get("v", np.zeros(J)).astype(np.float64).copy()
            dd = ws.get("d", np.zeros(I)).astype(np.float64).copy()
            theta_aging = np.asarray(
                ws.get("theta_aging", np.zeros(self.K_basis)), dtype=np.float64,
            ).copy()
            gamma = np.atleast_1d(
                np.asarray(ws.get("gamma", np.zeros(self.K_gamma)), dtype=np.float64),
            ).copy()
            if gamma.shape[0] != self.K_gamma:
                # wrong arity (scalar vs varying) — truncate / pad with zeros.
                tmp = np.zeros(self.K_gamma, dtype=np.float64)
                k = min(gamma.shape[0], self.K_gamma)
                tmp[:k] = gamma[:k]
                gamma = tmp
        else:
            if fcfg.init == "random":
                u = 0.01 * rng.standard_normal(I)
                v = 0.01 * rng.standard_normal(J)
            elif fcfg.init == "zeros":
                u = np.zeros(I); v = np.zeros(J)
            else:  # "mean" — u_i = mean response of athlete i, v_j = 0.
                usum = np.bincount(d.row_idx, weights=d.y, minlength=I)
                u = np.where(self.n_per_ath > 0,
                             usum / np.maximum(self.n_per_ath, 1.0), 0.0)
                v = np.zeros(J)
            dd = np.zeros(I)
            theta_aging = np.zeros(self.K_basis, dtype=np.float64)
            gamma = np.zeros(self.K_gamma, dtype=np.float64)

        if cfg.sigma2_init is not None:
            sigma2 = float(cfg.sigma2_init)
        else:
            r0 = d.y - u[d.row_idx] - v[d.col_idx]
            sigma2 = float(np.mean(r0 * r0)) + _EPS

        omega_d2 = float(cfg.omega_d2_init if cfg.omega_d2_init is not None
                         else cfg.omega_d2_fixed if cfg.omega_d2_fixed is not None
                         else _OMEGA_D2_INIT_DEFAULT)

        self.params = dict(
            u=u, v=v, d=dd,
            theta_aging=theta_aging, gamma=gamma,
            sigma2=sigma2, omega_d2=omega_d2,
            nu=float(cfg.nu),
        )

    # ------------------------------------------------------------------
    # Prediction (used by BaseModel.predict / residuals / log_lik)
    # ------------------------------------------------------------------
    def _aging_pred(self, theta_aging: np.ndarray) -> np.ndarray:
        # (K,) @ (K, N) -> (N,) via BLAS gemv; no (K, N) temporary.
        return theta_aging @ self.x_scalar[: self.K_basis]

    def _gamma_pred(self, gamma: np.ndarray) -> np.ndarray:
        return gamma @ self.x_scalar[self.K_basis :]

    def _yhat_terms(self) -> dict[str, np.ndarray]:
        d = self.data
        p = self.params
        out: dict[str, np.ndarray] = {}
        out["u"] = p["u"][d.row_idx]
        out["v"] = p["v"][d.col_idx]
        if self.config.use_phi12:
            out["aging"] = self._aging_pred(p["theta_aging"])
        if self.config.use_gamma:
            out["gamma"] = self._gamma_pred(p["gamma"])
        if self.config.use_d:
            out["d"] = p["d"][d.row_idx] * self.An_tilde_cell
        return out

    # ------------------------------------------------------------------
    # Objective: penalized Student-t log-likelihood
    # ------------------------------------------------------------------
    def _eb_penalty(
        self, d: np.ndarray | None = None, omega_d2: float | None = None,
    ) -> float:
        """EB log-prior term  -1/(2 omega_d^2) * sum_{i eligible} d_i^2.

        This is the penalty in the §7.1 objective. It is part of the function
        being maximized whenever the d_i block carries its N(0, omega_d^2) ridge
        — i.e. whenever use_d is on and any athlete is eligible — *regardless* of
        whether omega_d^2 is EB-estimated or frozen (the ridge is applied in the
        d update either way). Returns 0.0 when the d block is inactive.

        Note this is the quadratic prior term at the current omega_d^2, not the
        full type-II marginal (which would add the -1/2 * n_elig * log omega_d^2
        normalizer); see model_description_v2.md §7.1 / §7.9 for why the tracked
        objective uses the quadratic form.
        """
        if not (self.config.use_d and self.eligible_d.any()):
            return 0.0
        dd = self.params["d"] if d is None else d
        om2 = self.params["omega_d2"] if omega_d2 is None else omega_d2
        d_e = dd[self.eligible_d]
        return -0.5 * float(np.sum(d_e * d_e)) / max(om2, _EPS)

    # ------------------------------------------------------------------
    # APC gauge (§4.2 of model_description_v2.md)
    # ------------------------------------------------------------------
    def apply_apc_gauge_beta0(self) -> float:
        """Enforce the G1 (β=0) APC gauge in place. Returns c (the slope removed).

        The aging slope's linear-in-A^n component is not identified from the
        data alone: the transformation

            theta_aging ← theta_aging + delta(c)    [adds c*A^n to curve]
            v_j         ← v_j - c * t_j              [removes secular trend]
            u_i         ← u_i + c * f_i              [adds debut-cohort trend]

        leaves every cell prediction unchanged for any c, because
        A^n_ij = t_j - f_i exactly. G1 (β=0) fixes c so that the OLS slope of
        v_j on race date t_j is zero, pushing secular change into the aging
        curve.

        Safe to call at any point (per-iteration or post-fit). Re-imposes
        mean(v)=0 after adjusting v. Returns 0.0 and is a no-op when
        `use_phi12=False` (no aging curve to absorb the shift) or when the race
        dates span < _EPS fractional years (degenerate slice).
        """
        if not self.config.use_phi12:
            return 0.0

        v = self.params["v"]
        t = self._t_j   # (J,) race dates in fractional years

        # OLS slope c = cov(v, t) / var(t), centred.
        t_c = t - t.mean()
        denom = float(np.dot(t_c, t_c))
        if denom < _EPS:
            return 0.0
        c = float(np.dot(v, t_c) / denom)

        # v_j -= c * t_j  (removes secular trend)
        self.params["v"] = v - c * t

        # u_i += c * f_i  (adds debut-cohort trend)
        self.params["u"] = self.params["u"] + c * self._f_i

        # theta_aging: add c * A^n to the aging curve.
        # Each basis's first column is A^n (linear term); for orthogonalized poly
        # the coefficient-space delta is (c/scale)*R[:,0] — see §A.6.
        theta = self.params["theta_aging"].copy()
        cfg = self.config
        if cfg.basis_kind == "spline":
            # ncs_basis col 0 is literally A^n.
            theta[0] += c
        elif cfg.basis_kind == "poly":
            if cfg.orthogonalize:
                R, scale = self._basis_transform
                # delta s.t. B @ delta = c * A^n: derived from B = B_raw @ R^-1 * scale.
                theta += (c / scale) * R[:, 0]
            else:
                # raw monomials: col 0 = A^n.
                theta[0] += c
        self.params["theta_aging"] = theta

        # Re-impose mean(v) = 0 (level gauge, §4.1).
        shift = float(self.params["v"].mean())
        if abs(shift) > _EPS:
            self.params["v"] -= shift
            self.params["u"] += shift

        return c

    def log_lik(self) -> float:
        """Penalized objective: Student-t data loglik + EB log-prior on d_i.

        Overrides BaseModel.log_lik (data term only) so aic/bic/summary and
        FitResult.loglik_final all report the quantity the fitter actually
        maximizes. The held-out CV density (crossval.py) deliberately stays the
        data term only — a prior on parameters is not a likelihood of new data.
        """
        return super().log_lik() + self._eb_penalty()

    # ------------------------------------------------------------------
    # The fitter — fit() drives convergence; one iteration's block updates
    # live in _run_one_iteration() so solver variants (Anderson) reuse it.
    # ------------------------------------------------------------------
    def fit(self) -> FitResult:
        self._init_params()
        self._begin_fit()
        fcfg = self.fitter

        trace = self._new_trace()
        rss_prev = np.inf
        loglik_prev = -np.inf
        converged = False
        it = 0
        for it in range(1, fcfg.max_outer_iter + 1):
            pre = self._snapshot_blocks() if fcfg.record_trace else None
            rss = self._run_one_iteration()
            loglik = self._last_loglik
            if fcfg.record_trace:
                self._record_trace_row(trace, rss, pre)
            if fcfg.verbose and (it == 1 or it % 25 == 0):
                p = self.params
                print(f"  iter {it:4d}  rss={rss:.4e}  loglik={loglik:.4e}  "
                      f"sigma2={p['sigma2']:.3e}  om_d2={p['omega_d2']:.3e}")
            if self.iter_hook is not None:
                self.iter_hook(self, it, loglik, rss)
            if self._stop(loglik, loglik_prev, rss, rss_prev):
                converged = True
                break
            rss_prev = rss
            loglik_prev = loglik

        return self._finalize_fit(n_iter=it, converged=converged, trace=trace)

    def _stop(self, loglik: float, loglik_prev: float,
              rss: float, rss_prev: float) -> bool:
        tol = self.fitter.tol
        if self.fitter.stop_criterion == "loglik":
            return abs(loglik - loglik_prev) < tol * (abs(loglik) + _EPS)
        return abs(rss_prev - rss) < tol * (abs(rss) + _EPS)

    # ------------------------------------------------------------------
    # Trace helpers (shared with solver variants)
    # ------------------------------------------------------------------
    @staticmethod
    def _new_trace() -> dict[str, list[float]]:
        return {
            "rss": [], "loglik": [],
            "sigma2": [], "omega_d2": [],
            "max_du": [], "max_dv": [], "max_dd": [],
            "mean_du": [], "mean_dv": [], "mean_dd": [],
            "median_du": [], "median_dv": [], "median_dd": [],
        }

    def _snapshot_blocks(self) -> dict[str, np.ndarray]:
        p = self.params
        return {"u": p["u"].copy(), "v": p["v"].copy(), "d": p["d"].copy()}

    def _record_trace_row(
        self,
        trace: dict[str, list[float]],
        rss: float,
        pre: dict[str, np.ndarray] | None,
    ) -> None:
        p = self.params
        trace["rss"].append(rss)
        # Reuse the value cached in _run_one_iteration; do NOT inline a
        # self.log_lik() default (Python evaluates it unconditionally).
        cached = getattr(self, "_last_loglik", None)
        trace["loglik"].append(float(cached) if cached is not None else self.log_lik())
        trace["sigma2"].append(p["sigma2"])
        trace["omega_d2"].append(p["omega_d2"])
        if pre is None:
            for blk in ("u", "v", "d"):
                trace[f"max_d{blk}"].append(np.nan)
                trace[f"mean_d{blk}"].append(np.nan)
                trace[f"median_d{blk}"].append(np.nan)
        else:
            for blk in ("u", "v", "d"):
                diff = np.abs(p[blk] - pre[blk])
                trace[f"max_d{blk}"].append(float(diff.max()) if diff.size else np.nan)
                trace[f"mean_d{blk}"].append(float(diff.mean()) if diff.size else np.nan)
                trace[f"median_d{blk}"].append(float(np.median(diff)) if diff.size else np.nan)

    # ------------------------------------------------------------------
    # Shared fit machinery for solver variants
    # ------------------------------------------------------------------
    def _begin_fit(self) -> None:
        """Set per-fit scratch state. Call after _init_params()."""
        self._post_var = {"d": np.zeros(self.data.I)}
        self._irls_w = np.ones(self.data.N)

    def _finalize_fit(self, *, n_iter: int, converged: bool,
                       trace: dict[str, list[float]]) -> FitResult:
        self._residuals = self.data.y - self.predict()
        self.fit_result = FitResult(
            n_iter=n_iter,
            converged=converged,
            rss_final=float(np.sum(self._residuals ** 2)),
            loglik_final=self.log_lik(),
            trace=trace,
        )
        return self.fit_result

    # ------------------------------------------------------------------
    # The aging + gamma block solver
    # ------------------------------------------------------------------
    def _solve_aging_block(
        self, r: np.ndarray, w_rw: np.ndarray, theta: np.ndarray,
    ) -> None:
        """Update active entries of `theta` in place, patching `r` in place.

        `theta` is the full (K_blk,) vector laid out as in self.x_scalar.
          - "numba_cd": cd_inner_scalars (BCD within the block); poly only.
          - "direct":   one dense np.linalg.solve on the K_active Gram matrix.
        """
        if not self.scalar_use_mask.any():
            return

        # Fit the block with the (possibly capped) weight, but patch `r` over
        # every cell below so capped cells keep the extrapolated curve.
        w_blk = w_rw if self._aging_cell_w is None else w_rw * self._aging_cell_w

        if self.fitter.inner_solver == "numba_cd":
            cd_inner_scalars(
                r, w_blk, self.x_scalar,
                self.scalar_use_mask, theta,
                int(self.fitter.cd_inner_iters),
            )
            return
        if self.fitter.inner_solver != "direct":
            raise ValueError(f"unknown inner_solver={self.fitter.inner_solver!r}")

        mask = self.scalar_use_mask
        X = self._x_active                         # (K_act, N), cached in _precompute
        theta_old = theta[mask].copy()
        r_part = r + X.T @ theta_old               # partial residual
        WX = X * w_blk
        G = WX @ X.T                               # (K_act, K_act)
        b = WX @ r_part
        G.flat[:: G.shape[0] + 1] += _EPS          # tiny ridge for safety
        theta_new = np.linalg.solve(G, b)
        r -= X.T @ (theta_new - theta_old)
        theta[mask] = theta_new

    def _run_one_iteration(self) -> float:
        """One outer BCD iteration, in place on self.params. Returns RSS.

        Reads and writes self.params and self._post_var, so it is safe to call
        after an arbitrary state mutation (e.g. an Anderson mixing step).
        """
        d, cfg = self.data, self.config
        row, col = d.row_idx, d.col_idx
        I, J, N = d.I, d.J, d.N
        y = d.y

        p = self.params
        u, v, dd = p["u"], p["v"], p["d"]
        sigma2 = p["sigma2"]
        omega_d2 = p["omega_d2"]
        nu = p["nu"]
        theta = np.empty(self.K_blk, dtype=np.float64)
        theta[: self.K_basis] = p["theta_aging"]
        theta[self.K_basis :] = p["gamma"]

        omega_d2_locked = cfg.omega_d2_fixed is not None or cfg.freeze_eb_prior

        # --- residual refresh (kills float-drift) --------------------
        r = y - self._predict_full(u, v, dd, theta)

        # --- per-cell weights ---------------------------------------
        irls_w = irls_weights(r, sigma2, nu)
        if self.boot_w_cell is not None:
            irls_w = irls_w * self.boot_w_cell    # athlete bootstrap reweight
        w_base = irls_w / max(sigma2, _EPS)
        w_rw = w_base * self.w_ath_cell

        # --- v_j block (uses w_rw) ----------------------------------
        v_num = np.bincount(col, weights=w_rw * r, minlength=J)
        v_den = np.bincount(col, weights=w_rw, minlength=J)
        v_new = v + np.where(v_den > _EPS, v_num / np.maximum(v_den, _EPS), 0.0)
        r -= (v_new - v)[col]
        v = v_new

        # --- u_i block (uses w_base) --------------------------------
        u_num = np.bincount(row, weights=w_base * r, minlength=I)
        u_den = np.bincount(row, weights=w_base, minlength=I)
        u_new = u + np.where(u_den > _EPS, u_num / np.maximum(u_den, _EPS), 0.0)
        r -= (u_new - u)[row]
        u = u_new

        # --- {aging, gamma} block -----------------------------------
        self._solve_aging_block(r, w_rw, theta)

        # --- d_i block (uses w_base; EB ridge, eligible only) -------
        if cfg.use_d and self.eligible_d.any():
            Atc = self.An_tilde_cell
            wA = w_base * Atc
            d_num = np.bincount(row, weights=wA * r, minlength=I)
            d_den = np.bincount(row, weights=wA * Atc, minlength=I)
            ridge = 1.0 / max(omega_d2, _EPS)
            denom_d = d_den + ridge
            # Penalized Newton step. d_num = sum(w*Atc*r) is the *data* gradient
            # at the current dd (r already carries the -Atc*dd term); the ridge
            # adds the prior gradient -ridge*dd. Omitting -ridge*dd would make
            # the fixed point d_num=0 — the unpenalized MLE — with the ridge only
            # damping the step (see model_description_v2.md §7.4, d_i block).
            d_new = np.where(
                self.eligible_d,
                dd + (d_num - ridge * dd) / np.maximum(denom_d, _EPS),
                0.0,
            )
            r -= Atc * (d_new - dd)[row]
            dd = d_new
            self._post_var["d"] = np.where(
                self.eligible_d, 1.0 / np.maximum(denom_d, _EPS), 0.0,
            )

        # --- sigma^2 update (raw irls_w*r^2) ------------------------
        sigma2 = float(np.sum(irls_w * r * r) / max(N, 1)) + _EPS

        # --- EB update for omega_d2 ---------------------------------
        if cfg.use_d and self.eligible_d.any() and not omega_d2_locked:
            pv_d = self._post_var["d"][self.eligible_d]
            d_e = dd[self.eligible_d]
            new_om_d = float(np.mean(d_e * d_e + pv_d))
            if cfg.omega_d2_max is not None:
                new_om_d = min(new_om_d, float(cfg.omega_d2_max))
            omega_d2 = max(new_om_d, _EPS)
        elif cfg.omega_d2_fixed is not None:
            omega_d2 = float(cfg.omega_d2_fixed)

        # --- gauge: mean(v) = 0 (predict invariant under u+=c, v-=c) -
        if cfg.gauge == "mean_v_zero":
            c = float(v.mean())
            if abs(c) > _EPS:
                u += c
                v -= c

        self.params.update(
            u=u, v=v, d=dd,
            theta_aging=theta[: self.K_basis].copy(),
            gamma=theta[self.K_basis :].copy(),
            sigma2=sigma2, omega_d2=omega_d2,
        )
        self._irls_w = irls_w

        # APC gauge (§4.2): re-project onto the G1 manifold. Predictions are
        # invariant, so r and rss_val below are unaffected.
        if cfg.apc_gauge_beta0:
            self.apply_apc_gauge_beta0()

        # Penalized loglik from the in-loop residual (consistent with the updated
        # state), stashed for the stop check and the trace recorder. The EB prior
        # term must be included so the monitored objective matches §7.1 — without
        # it the loglik stop/trace track only the data fit, which a shrinkage step
        # legitimately worsens.
        rss_val = float(np.sum(r * r))
        data_ll = self._loglik_from_resid(r, sigma2, nu, N, rss=rss_val)
        self._last_loglik = data_ll + self._eb_penalty(dd, omega_d2)
        return rss_val

    @staticmethod
    def _loglik_from_resid(
        r: np.ndarray, sigma2: float, nu: float, N: int,
        *, rss: float | None = None,
    ) -> float:
        """Log-likelihood at the current state, from a precomputed residual."""
        if not np.isfinite(nu):
            if rss is None:
                rss = float(np.sum(r * r))
            return -0.5 * (
                N * np.log(2.0 * np.pi)
                + N * np.log(max(sigma2, 1e-300))
                + rss / max(sigma2, 1e-300)
            )
        c = N * (
            gammaln((nu + 1) / 2) - gammaln(nu / 2)
            - 0.5 * np.log(nu * np.pi * max(sigma2, 1e-300))
        )
        return float(c - 0.5 * (nu + 1)
                     * np.sum(np.log1p(r * r / (nu * max(sigma2, 1e-300)))))

    # Full predict from *passed-in* params, so the fit loop can build r before
    # writing params back. `theta` is the full (K_blk,) vector as in x_scalar.
    def _predict_full(self, u, v, dd, theta) -> np.ndarray:
        d = self.data
        out = u[d.row_idx] + v[d.col_idx]
        # Contiguous row-slices of x_scalar are views, so each `@` is a BLAS
        # gemv with only the (N,) output allocated — no (K, N) broadcast temp.
        if self.config.use_phi12 and self.K_basis > 0:
            out += theta[: self.K_basis] @ self.x_scalar[: self.K_basis]
        if self.config.use_gamma and self.K_gamma > 0:
            out += theta[self.K_basis :] @ self.x_scalar[self.K_basis :]
        if self.config.use_d:
            out += dd[d.row_idx] * self.An_tilde_cell
        return out

    # ------------------------------------------------------------------
    # Persistence of basis reconstruction artifacts
    # ------------------------------------------------------------------
    def _save_extra(self) -> dict[str, Any]:
        return {
            "basis_kind": self.config.basis_kind,
            "degree": self.config.degree,
            "orthogonalize": self.config.orthogonalize,
            "gamma_form": self.config.gamma_form,
            "K_basis": self.K_basis,
            "K_gamma": self.K_gamma,
            # data-derived, not recoverable from config alone:
            "poly_transform": self._basis_transform,   # (R, scale) or None
            "spline_knots": self._spline_knots,         # resolved knots or None
        }

    def _load_extra(self, extra: dict[str, Any]) -> None:
        if not extra:
            return
        self._basis_transform = extra.get("poly_transform")
        self._spline_knots = extra.get("spline_knots")

    # ------------------------------------------------------------------
    # AIC/BIC bookkeeping
    # ------------------------------------------------------------------
    def n_params_naive(self) -> int:
        I, J = self.data.I, self.data.J
        k = I + J - 1  # u + v under one gauge constraint
        if self.config.use_d:
            k += int(self.eligible_d.sum())
        if self.config.use_phi12:
            k += self.K_basis
        if self.config.use_gamma:
            k += self.K_gamma
        k += 1  # sigma2
        if self.config.use_d and not (
            self.config.omega_d2_fixed is not None or self.config.freeze_eb_prior
        ):
            k += 1
        return k

    def effective_dof(self) -> dict[str, float]:
        """Per-block effective d.o.f. via shrinkage.

        Non-EB blocks count their active parameters; the EB d_i block counts
        the trace of its hat matrix as sum(1 - shrinkage). The aging/gamma
        block has no EB prior, so its d.o.f. is its active column count.
        """
        I, J = self.data.I, self.data.J
        out: dict[str, float] = {
            "u_v": float(I + J - 1),                       # joint, one gauge
            "aging": float(self.K_basis if self.config.use_phi12 else 0),
            "gamma": float(self.K_gamma if self.config.use_gamma else 0),
        }
        if self.config.use_d and self.eligible_d.any() and "d" in self._post_var:
            om = max(self.params["omega_d2"], _EPS)
            ridge = 1.0 / om
            pv = np.maximum(self._post_var["d"][self.eligible_d], _EPS)
            signal = np.maximum(1.0 / pv - ridge, 0.0)
            out["d"] = float(np.sum(signal / (signal + ridge)))
        return out
