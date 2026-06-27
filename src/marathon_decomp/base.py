"""Contract every model variant satisfies.

Subclasses declare: NAME, PARAM_SPEC, CONFIG_CLS, _yhat_terms(), fit().
Base provides: predict, residuals, log_lik, aic, bic, save, load, summary.

Anything downstream (diagnostics, bootstrap, plotting) takes a BaseModel
and never touches model internals — so a new model variant gets all of
that for free once it implements the four required methods.
"""
from __future__ import annotations

import pickle
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

import numpy as np
from scipy.special import gammaln

from .data import FitData


@dataclass
class SaveSpec:
    params: bool = True
    residuals: bool = False
    irls_weights: bool = False
    posterior_var: bool = False
    convergence_trace: bool = False
    standard_errors: bool = False

    @classmethod
    def default(cls) -> "SaveSpec":
        return cls()

    @classmethod
    def full(cls) -> "SaveSpec":
        return cls(True, True, True, True, True, True)


@dataclass
class FitResult:
    """Outcome of a fit.

    `n_iter`, `converged`, `rss_final`, `loglik_final` are always populated
    and reflect the final state of the fit. `trace` contains per-iteration
    diagnostics (rss, loglik, sigma2, omega_*, per-block max/mean/median
    of |dx|) only when ``FitterConfig.record_trace=True``; otherwise it is
    an empty dict. Solver variants may add their own fields (e.g. Anderson
    populates ``anderson_accept_rate`` regardless of record_trace).

    Downstream code that consumes ``trace`` should defend against an empty
    dict.
    """

    n_iter: int
    converged: bool
    rss_final: float
    loglik_final: float
    trace: dict[str, list[float]] = field(default_factory=dict)


class BaseModel(ABC):
    """Abstract base for marathon decomposition models."""

    NAME: str = ""
    # Fitter version, persisted with every saved fit for provenance. Each
    # concrete model declares its own VERSION explicitly (rather than
    # inheriting) so variants can be bumped independently — e.g. a change to
    # the Anderson solver bumps ModelAnderson.VERSION without touching
    # Model.VERSION.
    VERSION: str = "0.0.0"
    PARAM_SPEC: dict[str, tuple] = {}
    CONFIG_CLS: type | None = None

    def __init__(self, data: FitData, config: Any):
        self.data = data
        self.config = config
        self.params: dict[str, Any] = {}
        self.fit_result: FitResult | None = None
        # populated by fit():
        self._residuals: np.ndarray | None = None
        self._irls_w: np.ndarray | None = None
        self._post_var: dict[str, np.ndarray] = {}
        # Optional per-iteration observer, fired once per outer iteration by
        # fit() as `iter_hook(self, it, loglik, rss)` after the iterate's state
        # and `_last_loglik` are resolved, before the stop check. None => no-op,
        # so it never alters fit output (no provenance/VERSION impact). Used to
        # snapshot per-iter factors (e.g. v) for convergence studies.
        self.iter_hook: Callable[["BaseModel", int, float, float], None] | None = None

    # --- subclass implements ---

    @abstractmethod
    def fit(self) -> FitResult: ...

    @abstractmethod
    def _yhat_terms(self) -> dict[str, np.ndarray]: ...

    def effective_dof(self) -> dict[str, float]:
        """Per-block effective d.o.f. Override for EB models."""
        return {}

    def n_params_naive(self) -> int:
        """Free-parameter count for naive AIC/BIC (no EB shrinkage)."""
        return 0

    def _save_extra(self) -> dict[str, Any]:
        """Model-specific data to persist for downstream reconstruction.

        Base returns nothing. Subclasses override to stash anything that is
        derived from data+config at fit time and not recoverable from the
        pickled `config`/`params` alone (e.g. the aging basis QR transform
        and resolved spline knots). Saved under payload['model_extra'].
        """
        return {}

    def _load_extra(self, extra: dict[str, Any]) -> None:
        """Restore whatever `_save_extra` persisted. No-op by default."""
        return None

    # --- base provides ---

    def predict(self) -> np.ndarray:
        terms = self._yhat_terms()
        if not terms:
            return np.zeros(self.data.N)
        out = np.zeros(self.data.N)
        for v in terms.values():
            out += v
        return out

    def residuals(self) -> np.ndarray:
        return self.data.y - self.predict()

    def log_lik(self) -> float:
        r = self.residuals()
        sigma2 = float(self.params.get("sigma2", np.var(r) + 1e-12))
        nu = float(self.params.get("nu", np.inf))
        N = self.data.N
        if not np.isfinite(nu):
            return -0.5 * (
                N * np.log(2.0 * np.pi)
                + N * np.log(max(sigma2, 1e-300))
                + float(np.sum(r * r)) / sigma2
            )
        c = N * (
            gammaln((nu + 1) / 2) - gammaln(nu / 2)
            - 0.5 * np.log(nu * np.pi * sigma2)
        )
        return float(c - 0.5 * (nu + 1) * np.sum(np.log1p(r * r / (nu * sigma2))))

    def aic(self, *, effective: bool = True) -> float:
        k = sum(self.effective_dof().values()) if effective else self.n_params_naive()
        return 2.0 * k - 2.0 * self.log_lik()

    def bic(self, *, effective: bool = True) -> float:
        k = sum(self.effective_dof().values()) if effective else self.n_params_naive()
        return k * np.log(self.data.N) - 2.0 * self.log_lik()

    def save(self, path: str | Path, *, what: SaveSpec | None = None) -> None:
        what = what or SaveSpec.default()
        path = Path(path)
        payload: dict[str, Any] = {
            "name": self.NAME,
            "fitter_version": self.VERSION,
            "config": self.config,
            "data_version": self.data.data_version,
            "spec": self.data.spec,
            "cache_key": self.data.cache_key,
            "fit_result": self.fit_result,
            "model_extra": self._save_extra(),
        }
        if what.params:
            payload["params"] = self.params
        if what.residuals:
            payload["residuals"] = self.residuals()
        if what.irls_weights and self._irls_w is not None:
            payload["irls_w"] = self._irls_w
        if what.posterior_var:
            payload["post_var"] = self._post_var
        if what.convergence_trace and self.fit_result is not None:
            payload["trace"] = self.fit_result.trace
        with path.open("wb") as f:
            pickle.dump(payload, f)

    @classmethod
    def load(cls, path: str | Path, data: FitData) -> "BaseModel":
        path = Path(path)
        with path.open("rb") as f:
            payload = pickle.load(f)
        m = cls(data, payload["config"])
        m.params = payload.get("params", {})
        m.fit_result = payload.get("fit_result")
        m._load_extra(payload.get("model_extra", {}))
        return m

    def summary(self) -> str:
        return (
            f"{self.NAME} v{self.VERSION}  I={self.data.I} J={self.data.J} N={self.data.N}  "
            f"loglik={self.log_lik():.2f}  aic={self.aic():.2f}  bic={self.bic():.2f}"
        )
