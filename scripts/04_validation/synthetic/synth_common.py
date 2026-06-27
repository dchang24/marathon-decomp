"""Self-contained synthetic-data validation library (scripts/04_validation/synthetic).

The experiment tests TWO claims that must be kept separate:

  1. NUMERICAL CONVERGENCE -- the solver reaches a stationary point of the
     objective J (docs/model_derivation.md S2.3). Diagnostics: the per-race
     normal-equation residual (S5.1) is ~0, the captured log-lik trace is
     monotone non-decreasing (EM theorem S4.1), the fit dominates the oracle
     (population) parameters, and ALS == Anderson at the fixed point.

  2. STATISTICAL RECOVERY -- that stationary point matches the known
     data-generating process (DGP). Diagnostics: gauge-aligned Pearson r and
     RMSE on u_i, v_j, d_i; scalar accuracy on sigma^2, the APC-invariant aging
     curvature (S7.2), and omega_d^2.

This module is deliberately STANDALONE. It does not import from `tests/`, so the
experiment directory is fully self-contained. The fast pytest gate in
tests/test_synthetic_recovery.py is left untouched; the DGP/metrics here are a
deliberate (documented) copy, free to diverge as the experiment grows.

Ground-truth magnitudes match the ALL_B_14-20 AxD nu6 reference fit (see
`GroundTruth`). All console output is ASCII-only (Windows cp1252).
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # scripts/

from marathon_decomp.aging import aging_curve_on_grid            # noqa: E402
from marathon_decomp.config import RESULTS_DIR                   # noqa: E402
from marathon_decomp.data import FitData, SliceSpec              # noqa: E402
from marathon_decomp.models.anderson import (                    # noqa: E402
    AndersonFitterConfig,
    ModelAnderson,
)
from marathon_decomp.models.model import (                       # noqa: E402
    FitterConfig,
    Model,
    ModelConfig,
)

_SEC_PER_YEAR = 365.25 * 86400.0   # same convention as data.py

OUT_ROOT = RESULTS_DIR / "validation" / "synthetic"


# ---------------------------------------------------------------------------
# Ground truth (ALL_B_14-20 AxD nu6 reference magnitudes)
# ---------------------------------------------------------------------------
# u: mean=9.6046, std=0.1957  (log-seconds; exp(9.6) ~ 14760 s ~ 4.1 h)
# v: std=0.0192               (race difficulty spread)
# d: std=0.0414               (per-year drift slope, eligible athletes only)
# sigma2: 0.00084             (sigma ~ 0.029 in log-time)
# omega_d2: 0.00192

@dataclass(frozen=True)
class GroundTruth:
    """Population latent-factor scales used to generate synthetic data."""

    u_mean: float = 9.60
    u_std: float = 0.20
    v_std: float = 0.020
    phi1: float = -0.010          # linear aging coef (raw poly-2 basis)
    phi2: float = 0.0022          # quadratic aging coef; curvature = 2*phi2
    gamma: float = 0.0015         # entry-age * career-age interaction
    d_std: float = float(np.sqrt(0.00192))   # ~0.044 per year
    sigma2: float = 0.00090       # sigma ~ 0.030
    omega_d2: float = 0.00192
    ae_lo: float = 18.0           # entry-age sampling band (years)
    ae_hi: float = 55.0

    @property
    def curvature(self) -> float:
        """d^2 f / dA_n^2 for the raw poly-2 basis (APC-gauge invariant)."""
        return 2.0 * self.phi2


GT = GroundTruth()


# ---------------------------------------------------------------------------
# Cell samplers -- two bipartite athlete-race designs
# ---------------------------------------------------------------------------

def sample_cells_random(
    I: int, J: int, target_N: int, rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Brinker's missingness model: a CONNECTED random bipartite graph.

    Guarantees every athlete has >=1 cell and every race has >=1 athlete (one
    forced valid entry per row and per column, exactly as in Brinker 2025
    sec. 4.2), then fills the rest at random to hit target_N. Degrees are
    binomial -> heterogeneous field sizes, like the real data.
    """
    cells: set[tuple[int, int]] = set()
    pa = rng.permutation(I)
    pj = rng.permutation(J)
    for i in range(I):
        cells.add((int(pa[i]), int(pj[i % J])))
    for j in range(J):
        cells.add((int(pa[j % I]), int(pj[j])))
    while len(cells) < target_N:
        cells.add((int(rng.integers(0, I)), int(rng.integers(0, J))))
    arr = np.array(sorted(cells), dtype=np.int64)
    return arr[:, 0], arr[:, 1]


def _is_connected(row: np.ndarray, col: np.ndarray, I: int, J: int) -> bool:
    """Union-find: is the athlete-race bipartite graph a single component?

    Factors are only jointly identified within a connected component (S7.1), so
    a disconnected design silently destroys cross-block v_j comparability.
    Athletes are nodes 0..I-1, races I..I+J-1.
    """
    parent = np.arange(I + J)

    def find(x: int) -> int:
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:
            parent[x], x = root, parent[x]
        return root

    for a, b in zip(row.tolist(), (col + I).tolist()):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    return len({find(x) for x in range(I + J)}) == 1


def sample_cells_staggered(
    I: int, J: int, k: int, rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Staggered-entry design: each athlete debuts at a random race and runs k
    races from there onward, so DEBUT DATE VARIES across athletes.

    Required to identify the aging curve: aging(A_n)=aging(t_j - b_i) is only
    separable from a per-race effect when the debut b_i varies across athletes
    (the age/cohort axis of the APC problem, S7.2). The balanced/random designs
    collapse debut to race 0 at high fill (everyone's first race is the first
    race), leaving aging fully confounded with v_j. Entry e_i ~ U[0, J-k] keeps
    row degree exactly k; later races draw from a larger eligible pool, so field
    sizes have a mild upward gradient (the realism cost of cohort variation).
    """
    k = int(max(1, min(k, J)))
    rows = np.repeat(np.arange(I, dtype=np.int64), k)
    cols = np.empty(I * k, dtype=np.int64)
    hi = max(1, J - k + 1)
    for _ in range(8):
        ok = True
        for i in range(I):
            e = int(rng.integers(0, hi))
            cols[i * k:(i + 1) * k] = e + rng.permutation(J - e)[:k]
        if _is_connected(rows, cols, I, J):
            ok = True
            break
        ok = False
    if not ok:
        raise RuntimeError(f"staggered sampler disconnected at I={I},J={J},k={k}")
    return rows, cols


def sample_cells_balanced(
    I: int, J: int, k: int, rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Symmetric "equal-scale" design: every athlete runs EXACTLY k distinct
    races (chosen uniformly at random), so every race draws ~ I*k/J runners.

    Vectorized: per-athlete, take the k smallest of J random keys -> k distinct
    uniform races. Row degree is exactly k; column degrees are binomial with
    mean I*k/J and sd ~ sqrt(I*k/J*(1-k/J)) (e.g. 200 +- 10 at I=400,J=40,k=20)
    -- "roughly similar" field sizes, the requested simplification of Brinker's
    heterogeneous design, while staying connected and well-mixed (no block
    structure). Connectivity is verified; a rare disconnect retries with fresh
    keys before raising.
    """
    k = int(max(1, min(k, J)))
    rows = np.repeat(np.arange(I, dtype=np.int64), k)
    for _ in range(8):
        keys = rng.random((I, J))
        cols = np.argsort(keys, axis=1)[:, :k].reshape(-1).astype(np.int64)
        if _is_connected(rows, cols, I, J):
            return rows, cols
    raise RuntimeError(
        f"balanced sampler could not produce a connected graph at "
        f"I={I}, J={J}, k={k} (k too small relative to J?)."
    )


# ---------------------------------------------------------------------------
# Synthetic data factory (faithful forward pass of model_derivation eq. 2)
# ---------------------------------------------------------------------------

@dataclass
class SynthTruth:
    """The realized ground-truth arrays for one synthetic dataset."""

    u: np.ndarray
    v: np.ndarray
    d: np.ndarray
    signal: np.ndarray
    terms: str
    gt: GroundTruth = field(default_factory=lambda: GT)


def make_synthetic(
    I: int = 400,
    J: int = 40,
    fill: float = 0.5,
    *,
    sampler: str = "balanced",        # "balanced" | "random"
    terms: str = "full",              # "full" | "rank1"
    noise_nu: float = float("inf"),
    seed: int = 0,
    gt: GroundTruth = GT,
) -> tuple[FitData, SynthTruth]:
    """Build a synthetic FitData from the model's own forward pass.

    `terms="full"` includes aging (phi1,phi2), the entry-age gamma fan, and the
    per-athlete drift d_i. `terms="rank1"` is the bare baseline u_i + v_j +
    noise (Brinker's rank-1 NMF analogue), used for the clean sparsity sweep.

    Noise is drawn as the exact Student-t scale mixture (a Gaussian whose
    precision is Gamma(nu/2, nu/2); model_derivation S3), so noise_nu finite
    yields genuine t_nu(0, sigma^2) residuals; noise_nu=inf is Gaussian.

    Debut date per athlete = min race date over that athlete's cells (the
    FitData convention), so A_n >= 0 everywhere without filtering.
    """
    rng = np.random.default_rng(seed)
    full = terms == "full"

    race_dates = pd.date_range("2015-01-01", "2024-12-31", periods=J).to_numpy()
    t_j = race_dates.astype("int64") * (1e-9 / _SEC_PER_YEAR)

    if sampler == "balanced":
        row, col = sample_cells_balanced(I, J, int(round(fill * J)), rng)
    elif sampler == "staggered":
        row, col = sample_cells_staggered(I, J, int(round(fill * J)), rng)
    elif sampler == "random":
        target_N = max(int(round(fill * I * J)), I + J + 10)
        row, col = sample_cells_random(I, J, target_N, rng)
    else:
        raise ValueError(f"unknown sampler={sampler!r}")
    N = len(row)

    f_i = np.full(I, np.inf)
    np.minimum.at(f_i, row, t_j[col])
    A_n = t_j[col] - f_i[row]

    A_e = rng.uniform(gt.ae_lo, gt.ae_hi, I)
    Ae_c_cell = (A_e - float(A_e.mean()))[row]

    n_per = np.bincount(row, minlength=I).astype(np.float64)
    An_bar_i = np.where(
        n_per > 0,
        np.bincount(row, weights=A_n, minlength=I) / np.maximum(n_per, 1.0),
        0.0,
    )
    An_tilde = A_n - An_bar_i[row]

    u_true = rng.normal(gt.u_mean, gt.u_std, I)
    v_raw = rng.normal(0.0, gt.v_std, J)
    v_true = v_raw - v_raw.mean()
    d_true = rng.normal(0.0, gt.d_std, I) if full else np.zeros(I)

    signal = u_true[row] + v_true[col]
    if full:
        signal = (
            signal
            + gt.phi1 * A_n + gt.phi2 * A_n ** 2
            + gt.gamma * Ae_c_cell * A_n
            + d_true[row] * An_tilde
        )

    sigma = float(np.sqrt(gt.sigma2))
    if np.isfinite(noise_nu):
        tau = rng.gamma(noise_nu / 2.0, 2.0 / noise_nu, N)
        noise = rng.normal(0.0, sigma / np.sqrt(tau))
    else:
        noise = rng.normal(0.0, sigma, N)
    y = signal + noise

    fd = FitData(
        row_idx=row.astype(np.int64), col_idx=col.astype(np.int64),
        y=y.astype(np.float64), A_n=A_n.astype(np.float64),
        A_e=A_e.astype(np.float64),
        athlete_ids=np.arange(I, dtype=np.int64),
        athlete_sex=np.array(["M"] * I, dtype=object),
        athlete_country=np.array(["GB"] * I, dtype=object),
        race_ids=np.arange(J, dtype=np.int64),
        race_series=np.array(["synth"] * J, dtype=object),
        race_country=np.array(["GB"] * J, dtype=object),
        race_date=race_dates,
        I=I, J=J, N=N,
        spec=SliceSpec(sex="M"),
        data_version="synthetic", cache_key="synthetic",
        response_kind="log_time",
    )
    truth = SynthTruth(u=u_true, v=v_true, d=d_true, signal=signal, terms=terms, gt=gt)
    return fd, truth


# ---------------------------------------------------------------------------
# Model / fitter config helpers
# ---------------------------------------------------------------------------

def model_config(terms: str = "full", nu: float = float("inf"), **overrides) -> ModelConfig:
    """Model config for the synthetic study.

    `terms="rank1"` disables the aging/gamma/drift blocks so the fit matches the
    rank-1 DGP. apc_gauge_beta0=False throughout: the APC gauge is tested
    separately, and per-iteration gauge shifts perturb the convergence trace.
    """
    full = terms == "full"
    cfg = ModelConfig(
        use_phi12=full, use_gamma=full, use_d=full,
        basis_kind="poly", degree=2, orthogonalize=True,
        gamma_form="scalar",
        omega_d2_max=0.05,
        apc_gauge_beta0=False,
        runner_weight_n0=1.0,
        nu=nu,
    )
    for key, val in overrides.items():
        setattr(cfg, key, val)
    return cfg


def anderson_fitter(**overrides) -> AndersonFitterConfig:
    cfg = AndersonFitterConfig(max_outer_iter=500, tol=1e-7, verbose=0)
    for key, val in overrides.items():
        setattr(cfg, key, val)
    return cfg


def als_fitter(**overrides) -> FitterConfig:
    cfg = FitterConfig(max_outer_iter=1000, tol=1e-7, verbose=0)
    for key, val in overrides.items():
        setattr(cfg, key, val)
    return cfg


def fit_with_trace(model) -> list[float]:
    """Attach an iter_hook that records the log-lik per outer iteration, fit,
    and return the trace. Used for the monotone-objective convergence check.
    """
    trace: list[float] = []
    model.iter_hook = lambda m, it, ll, rss: trace.append(float(ll))
    model.fit()
    return trace


# ---------------------------------------------------------------------------
# Recovery metrics (gauge-aware; model_derivation S7)
# ---------------------------------------------------------------------------

def _pearson(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 2 or np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def _poly_resid(x: np.ndarray, *bases_degrees) -> np.ndarray:
    """Residual of x after OLS on [1] + {base**1..base**deg for each (base, deg)}.

    Used to project v/u onto the identified subspace. For a degree-K aging
    basis, v_j is identified only up to a degree-K polynomial in race date
    (the APC gauge S7.2 generalized: aging(A_n)=aging(t_j - b_i) injects a
    degree-K polynomial in t_j into the per-race direction), and u_i symmetric
    in debut date. Removing that polynomial from BOTH fit and truth compares
    them in the same gauge.
    """
    cols = [np.ones_like(x, dtype=np.float64)]
    for base, deg in bases_degrees:
        for p in range(1, deg + 1):
            cols.append(np.asarray(base, dtype=np.float64) ** p)
    X = np.column_stack(cols)
    beta = np.linalg.lstsq(X, x.astype(np.float64), rcond=None)[0]
    return x - X @ beta


def _race_year_centered(fd: FitData) -> np.ndarray:
    tj = fd.race_date.astype("int64") * (1e-9 / _SEC_PER_YEAR)
    return tj - tj.mean()


def _athlete_debut_mean_centered(fd: FitData) -> tuple[np.ndarray, np.ndarray]:
    """Per-athlete debut (min) and mean race-year, both centered."""
    tj = fd.race_date.astype("int64") * (1e-9 / _SEC_PER_YEAR)
    debut = np.full(fd.I, np.inf)
    np.minimum.at(debut, fd.row_idx, tj[fd.col_idx])
    n = np.bincount(fd.row_idx, minlength=fd.I).astype(float)
    mean_t = (np.bincount(fd.row_idx, weights=tj[fd.col_idx], minlength=fd.I)
              / np.maximum(n, 1.0))
    return debut - debut.mean(), mean_t - mean_t.mean()


def _gaussian_loglik(y, yhat, sigma2) -> float:
    r = y - yhat
    N = len(y)
    return float(-0.5 * (N * np.log(2 * np.pi * sigma2) + np.sum(r * r) / sigma2))


def aging_curvature(model) -> float:
    """Central-difference d^2 f/dA_n^2 at A_n=2 yr (cancels the APC-free linear
    term; model_derivation S7.2). True value = gt.curvature = 2*phi2.
    """
    grid = np.array([1.90, 2.00, 2.10], dtype=np.float64)
    f = aging_curve_on_grid(model, grid)
    return float((f[2] - 2.0 * f[1] + f[0]) / 0.01)


def recovery_metrics(model, fd: FitData, truth: SynthTruth) -> dict:
    """Gauge-aligned recovery of every identified quantity.

    The headline r_v/r_u are GAUGE-FIXED (the only honest recovery target):
      v: residualized on a degree-K polynomial in race date (K = aging basis
         degree; 0 when aging is off) -- the APC gauge S7.2 generalized.
      u: residualized on a degree-K polynomial in debut date + mean race date
         (the cohort/drift tilts S7.2-S7.3).
    Raw (centre-only) r_v_raw/r_u_raw are kept so the gauge gap is explicit.
    d: eligible set only, centered (only differences are identified, S7.3).
    aging: APC-invariant curvature (S7.2). Plus Brinker MRE on exp(y).
    """
    full = truth.terms == "full"
    cfg = model.config
    K = cfg.degree if (full and getattr(cfg, "use_phi12", False)
                       and cfg.basis_kind == "poly") else 0

    u_fit = np.asarray(model.params["u"], dtype=np.float64)
    v_fit = np.asarray(model.params["v"], dtype=np.float64)

    yr = _race_year_centered(fd)
    debut_c, meant_c = _athlete_debut_mean_centered(fd)

    # gauge-fixed (headline)
    vf_g, vt_g = _poly_resid(v_fit, (yr, K)), _poly_resid(truth.v, (yr, K))
    uf_g = _poly_resid(u_fit, (debut_c, K), (meant_c, 1 if full else 0))
    ut_g = _poly_resid(truth.u, (debut_c, K), (meant_c, 1 if full else 0))

    yhat = model.predict()
    mre = float(np.mean(np.abs(np.expm1(yhat - fd.y))))   # |exp(yhat)/exp(y) - 1|

    out = {
        "r_v": _pearson(vt_g, vf_g),
        "rmse_v": float(np.sqrt(np.mean((vt_g - vf_g) ** 2))),
        "r_u": _pearson(ut_g, uf_g),
        "rmse_u": float(np.sqrt(np.mean((ut_g - uf_g) ** 2))),
        "r_v_raw": _pearson(truth.v - truth.v.mean(), v_fit - v_fit.mean()),
        "r_u_raw": _pearson(truth.u - truth.u.mean(), u_fit - u_fit.mean()),
        "sigma2_ratio": float(model.params["sigma2"] / truth.gt.sigma2),
        "mre": mre,
    }

    if full:
        elig = model.eligible_d
        d_fit = np.asarray(model.params["d"], dtype=np.float64)
        if elig.any():
            d_t = truth.d[elig] - truth.d[elig].mean()
            d_f = d_fit[elig] - d_fit[elig].mean()
            out["r_d"] = _pearson(d_t, d_f)
            out["rmse_d"] = float(np.sqrt(np.mean((d_t - d_f) ** 2)))
            out["d_shrink_slope"] = (
                float(np.polyfit(truth.d[elig], d_fit[elig], 1)[0])
                if np.std(truth.d[elig]) > 0 else float("nan")
            )
        else:
            out["r_d"] = out["rmse_d"] = out["d_shrink_slope"] = float("nan")
        out["curvature_err"] = aging_curvature(model) - truth.gt.curvature
        out["omega_d2_ratio"] = float(model.params["omega_d2"] / truth.gt.omega_d2)

    return out


# ---------------------------------------------------------------------------
# Convergence diagnostics (objective stationarity; model_derivation S2-S5)
# ---------------------------------------------------------------------------

def convergence_diagnostics(model, fd: FitData, truth: SynthTruth,
                            trace: list[float] | None = None) -> dict:
    """Numerical-convergence checks, independent of the ground truth value.

    stationarity_v : max_j |sum_i w_ij r_ij| / sum_i w_ij -- the size of the
                     v-block update step at the fixed point (S5.1 normal eqn).
                     A converged BCD has this ~ 0.
    oracle_margin  : log_lik(fit) - log_lik(oracle DGP params). MLE on a noisy
                     sample must beat the population-optimal oracle, so > 0.
    mono_min_step  : min successive delta of the log-lik trace; >= ~0 confirms
                     the EM monotonicity theorem (S4.1) held numerically.
    """
    res = model.fit_result
    r = model.residuals()
    nu = float(model.params.get("nu", np.inf))
    sigma2 = float(model.params["sigma2"])
    w = (nu + 1.0) / (nu + r * r / sigma2) if np.isfinite(nu) else np.ones(fd.N)

    num = np.bincount(fd.col_idx, weights=w * r, minlength=fd.J)
    den = np.bincount(fd.col_idx, weights=w, minlength=fd.J)
    stationarity_v = float(np.max(np.abs(num) / np.maximum(den, 1e-300)))

    if np.isfinite(nu):
        from scipy.special import gammaln
        rr = fd.y - truth.signal
        c = fd.N * (gammaln((nu + 1) / 2) - gammaln(nu / 2)
                    - 0.5 * np.log(nu * np.pi * truth.gt.sigma2))
        ll_oracle = float(c - 0.5 * (nu + 1) * np.sum(
            np.log1p(rr * rr / (nu * truth.gt.sigma2))))
    else:
        ll_oracle = _gaussian_loglik(fd.y, truth.signal, truth.gt.sigma2)

    out = {
        "converged": bool(res.converged),
        "n_iter": int(res.n_iter),
        "stationarity_v": stationarity_v,
        "oracle_margin": float(model.log_lik() - ll_oracle),
    }
    if trace is not None and len(trace) >= 2:
        steps = np.diff(np.asarray(trace, dtype=np.float64))
        out["mono_min_step"] = float(steps.min())
    return out


# ---------------------------------------------------------------------------
# Reporting + table I/O
# ---------------------------------------------------------------------------

class Report:
    """Console + markdown sink (ASCII only). Same house style as delta_common."""

    def __init__(self, title: str):
        self.lines: list[str] = []
        self.h1(title)

    def _emit(self, s: str = "") -> None:
        print(s)
        self.lines.append(s)

    def h1(self, s: str) -> None:
        self._emit(); self._emit(f"# {s}"); self._emit()

    def h2(self, s: str) -> None:
        self._emit(); self._emit(f"## {s}"); self._emit()

    def line(self, s: str = "") -> None:
        self._emit(s)

    def table(self, df: pd.DataFrame, floatfmt: str = "{:.4g}") -> None:
        """Render `df` to the console (to_string) and as a markdown pipe-table.

        Dependency-free (no `tabulate`): NaN/inf shown as '-'. Same house style
        as scripts/04_validation/d_i/delta_common.Report.
        """
        def fmt(x):
            if isinstance(x, (float, np.floating)):
                return "-" if not np.isfinite(x) else floatfmt.format(x)
            return str(x)

        disp = df.copy()
        for c in disp.columns:
            disp[c] = disp[c].map(fmt)
        with pd.option_context("display.width", 200, "display.max_columns", 60):
            print(disp.to_string(index=False), flush=True)
        cols = list(disp.columns)
        self.lines.append("| " + " | ".join(cols) + " |")
        self.lines.append("| " + " | ".join("---" for _ in cols) + " |")
        for _, r in disp.iterrows():
            self.lines.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
        self.lines.append("")

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(self.lines) + "\n", encoding="utf-8")


def upsert_csv(df_new: pd.DataFrame, path: Path, key: list[str]) -> pd.DataFrame:
    """Merge df_new into the CSV at `path`, replacing rows matching on `key`."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        old = pd.read_csv(path)
        merged = pd.concat([old, df_new], ignore_index=True)
        merged = merged.drop_duplicates(subset=key, keep="last").reset_index(drop=True)
    else:
        merged = df_new
    merged.to_csv(path, index=False)
    return merged
