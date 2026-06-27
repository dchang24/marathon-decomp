"""ModelAnderson — the same model as Model, with Anderson acceleration
wrapping the outer fixed-point iteration.

The outer loop is a fixed-point iteration  Theta_{k+1} = G(Theta_k) where G is
one full BCD sweep (Model._run_one_iteration). BCD fixed points converge
linearly; Anderson replaces the plain step with an extrapolation over the last
m iterates, turning linear convergence into roughly q-superlinear.

Per iteration:
    x_k = pack(state)
    apply G:  state <- G(x_k);  g_k = pack(state)
    push (x_k, g_k) to history (max length m)
    if history >= 2:
        solve   min ||F alpha||^2   s.t.   sum alpha = 1     (F = G - X residuals)
        propose x_prop = G_mat @ alpha
        if loglik(x_prop) >= loglik_plain - slack*|loglik_plain|:  keep proposal
        else:                                                      revert to g_k

The safeguard gates on the actual objective (the Student-t log-likelihood),
not the unweighted RSS: under finite nu an RSS gate rewards pulling outliers
in, which stalls the solver short of the optimum. In the Gaussian limit the
loglik gate reduces to the RSS one.

State packing: all parameters concatenated into one vector; positive scalars
(sigma^2, omega_d^2) live in log-space so the mix can't propose negatives.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..base import FitResult
from .model import FitterConfig, Model, ModelConfig

_EPS = 1e-12


@dataclass
class AndersonFitterConfig(FitterConfig):
    """FitterConfig + Anderson-specific knobs."""

    anderson_m: int = 5                       # history length
    anderson_start_iter: int = 3              # plain BCD before mixing
    anderson_safeguard_slack: float = 0.0     # accept if loglik_prop >= loglik_plain - slack*|loglik_plain|
    anderson_reset_on_reject: bool = True
    anderson_reg: float = 1e-10               # ridge on normal-eqn matrix
    # relaxation: x_prop = (1-beta)*sum a_i x_i + beta*sum a_i g_i. beta=1 is
    # full (undamped) Anderson; beta in (0,1) damps toward the plain BCD average.
    anderson_beta: float = 1.0


class ModelAnderson(Model):
    NAME = "model_anderson"
    # 1.3.0: inherits the _t_j units fix that re-enables the beta=0 APC gauge
    # (Model._setup_aging). Fits saved by <= 1.2.0 are off-manifold.
    VERSION = "1.3.0"

    def __init__(
        self,
        data,
        config: ModelConfig | None = None,
        fitter: AndersonFitterConfig | None = None,
    ):
        super().__init__(data, config, fitter or AndersonFitterConfig())

    # ------------------------------------------------------------------
    # Pack / unpack
    # ------------------------------------------------------------------
    def _pack(self) -> np.ndarray:
        p = self.params
        return np.concatenate([
            p["u"], p["v"], p["d"],
            np.asarray(p["theta_aging"], dtype=np.float64),
            np.atleast_1d(np.asarray(p["gamma"], dtype=np.float64)),
            np.array([
                np.log(max(p["sigma2"], _EPS)),
                np.log(max(p["omega_d2"], _EPS)),
            ]),
        ])

    def _unpack(self, vec: np.ndarray) -> None:
        I, J = self.data.I, self.data.J
        Kb, Kg = self.K_basis, self.K_gamma
        p = self.params
        i = 0
        p["u"] = vec[i:i + I].copy(); i += I
        p["v"] = vec[i:i + J].copy(); i += J
        p["d"] = vec[i:i + I].copy(); i += I
        p["theta_aging"] = vec[i:i + Kb].copy(); i += Kb
        p["gamma"] = vec[i:i + Kg].copy(); i += Kg
        p["sigma2"] = float(np.exp(vec[i])); i += 1
        p["omega_d2"] = float(np.exp(vec[i])); i += 1

    def _rss_at_state(self) -> float:
        """RSS at the current params, also caching `_last_loglik` so callers
        that change state via mixing don't need a separate predict()."""
        r = self.data.y - self.predict()
        rss = float(np.sum(r * r))
        p = self.params
        # Penalized objective so the safeguard gates on the §7.1 function: a
        # proposal that shrinks d_i trades data fit for prior fit, and a
        # data-only gate would wrongly reject it.
        data_ll = self._loglik_from_resid(
            r, p["sigma2"], p["nu"], self.data.N, rss=rss,
        )
        self._last_loglik = data_ll + self._eb_penalty()
        return rss

    # ------------------------------------------------------------------
    # Anderson outer loop
    # ------------------------------------------------------------------
    def fit(self) -> FitResult:
        self._init_params()
        self._begin_fit()
        fcfg: AndersonFitterConfig = self.fitter

        X_hist: list[np.ndarray] = []
        G_hist: list[np.ndarray] = []
        trace = self._new_trace()
        trace["anderson_accepted"] = []

        rss_prev = np.inf
        loglik_prev = -np.inf
        converged = False
        n_accept = 0
        it = 0
        for it in range(1, fcfg.max_outer_iter + 1):
            pre = self._snapshot_blocks() if fcfg.record_trace else None
            x_k = self._pack()
            rss_plain = self._run_one_iteration()
            loglik_plain = float(self._last_loglik)  # cached at post-G state
            g_k = self._pack()

            X_hist.append(x_k); G_hist.append(g_k)
            if len(X_hist) > fcfg.anderson_m:
                X_hist.pop(0); G_hist.pop(0)

            rss_chosen = rss_plain
            accepted = False
            if it >= fcfg.anderson_start_iter and len(X_hist) >= 2:
                accepted, rss_chosen = self._try_anderson_mix(
                    X_hist, G_hist, g_k, rss_plain, loglik_plain, fcfg,
                )
                if accepted:
                    n_accept += 1
                elif fcfg.anderson_reset_on_reject:
                    X_hist.clear(); G_hist.clear()

            # Resolve `_last_loglik` for the chosen state before recording the
            # trace: on accept it reflects the proposal (set by _rss_at_state);
            # on reject the proposal value is stale, so use the plain-G value.
            if accepted:
                loglik_chosen = float(self._last_loglik)
            else:
                loglik_chosen = loglik_plain
                self._last_loglik = loglik_chosen

            if fcfg.record_trace:
                self._record_trace_row(trace, rss_chosen, pre)
                trace["anderson_accepted"].append(accepted)

            if fcfg.verbose and (it == 1 or it % 25 == 0):
                tag = "A" if accepted else "."
                p = self.params
                print(f"  iter {it:4d} [{tag}]  rss={rss_chosen:.4e}  "
                      f"loglik={loglik_chosen:.4e}  "
                      f"sigma2={p['sigma2']:.3e}  om_d2={p['omega_d2']:.3e}")

            if self.iter_hook is not None:
                self.iter_hook(self, it, loglik_chosen, rss_chosen)

            if self._stop(loglik_chosen, loglik_prev, rss_chosen, rss_prev):
                converged = True
                break
            rss_prev = rss_chosen
            loglik_prev = loglik_chosen

        result = self._finalize_fit(n_iter=it, converged=converged, trace=trace)
        result.trace.setdefault("anderson_accept_rate", []).append(
            n_accept / max(it, 1)
        )
        return result

    # ------------------------------------------------------------------
    def _try_anderson_mix(
        self,
        X_hist: list[np.ndarray],
        G_hist: list[np.ndarray],
        g_k: np.ndarray,
        rss_plain: float,
        loglik_plain: float,
        fcfg: AndersonFitterConfig,
    ) -> tuple[bool, float]:
        """Attempt one Anderson mixing step. Returns (accepted, rss_chosen).

        On accept: state is left at the mixed proposal.
        On reject: state is restored to g_k (the plain-G post-step state).
        """
        m = len(X_hist)
        # Columns: most-recent first.
        F = np.column_stack([G_hist[m - 1 - i] - X_hist[m - 1 - i]
                             for i in range(m)])
        G_mat = np.column_stack([G_hist[m - 1 - i] for i in range(m)])

        # Reparameterise: alpha_0 = 1 - sum(beta), alpha_{1:} = beta.
        # min ||f_0 + dF @ beta||^2  =>  beta = -(dF^T dF)^{-1} dF^T f_0.
        f0 = F[:, 0]
        dF = F[:, 1:] - f0[:, None]
        A = dF.T @ dF
        diag_scale = np.maximum(np.diag(A).max(), 1.0) if A.size else 1.0
        A_reg = A + fcfg.anderson_reg * diag_scale * np.eye(A.shape[0])
        rhs = -dF.T @ f0
        try:
            beta = np.linalg.solve(A_reg, rhs)
        except np.linalg.LinAlgError:
            return False, rss_plain    # state still at g_k

        alpha = np.concatenate([[1.0 - beta.sum()], beta])
        b = fcfg.anderson_beta
        if b == 1.0:
            x_prop = G_mat @ alpha
        else:
            X_mat = np.column_stack([X_hist[m - 1 - i] for i in range(m)])
            x_prop = X_mat @ alpha + b * (G_mat - X_mat) @ alpha

        # Safeguard on the (Student-t) log-likelihood — see module docstring.
        self._unpack(x_prop)
        rss_prop = self._rss_at_state()           # also caches _last_loglik
        loglik_prop = float(self._last_loglik)
        if loglik_prop >= loglik_plain - fcfg.anderson_safeguard_slack * abs(loglik_plain):
            return True, rss_prop
        self._unpack(g_k)                          # reject: restore plain-G state
        return False, rss_plain
