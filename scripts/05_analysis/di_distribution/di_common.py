"""Shared loader + descriptors for the d_i (career-drift) distribution study.

Runs off the registered production **AxD** point fit (`full_nu8p00_best`) -- no
refit, no bootstrap file needed. Per *eligible* athlete (n_i >= 3) we read the EB
drift slope `d_i` (`params["d"]`) and rebuild its EB posterior variance
`post_var_i` from residuals exactly as the drift diagnostics do (`post_var` is not
persisted): `s_den_i = sum_cells w_base * ac^2`, `w_base = irls_w / sigma2`,
`ac = A_n - mean_i A_n`, `post_var_i = 1 / (s_den_i + 1/omega_d2)`.

Sign convention (project-wide): `d_i < 0` = the athlete runs *faster* than the
shared population aging trajectory as career age grows -> an **improver**;
`d_i > 0` = a **decliner**. Likewise `u_i < 0` = a faster athlete.

IDENTIFICATION CAVEAT (see the dir README + [[aging-slope-apc-nonidentification]]):
adding a constant c to every d_i is absorbed by the aging block's linear
coefficient (+c*A_n) and u_i (-c*mean_i A_n). So the **level/mean of d_i is not
identified by the likelihood** -- it is pinned only by the EB prior (centers it at
~0). Therefore `frac_improver` and `corr(d, u)` are **gauge-/prior-anchored**;
the **disattenuated SD and skew are gauge-invariant** (a location shift changes
neither) and are the descriptors to lean on.
"""
from __future__ import annotations

import functools
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # scripts/

from marathon_decomp import load_slice, registry                # noqa: E402
from marathon_decomp.config import RESULTS_DIR                  # noqa: E402
from marathon_decomp.kernels import irls_weights               # noqa: E402

from baseline_common import slices as S                          # noqa: E402
from drift_common.fitting import drift_cfg, drift_stem           # noqa: E402

OUT_ROOT = RESULTS_DIR / "analysis" / "di_distribution"
MODELS_ROOT = RESULTS_DIR / "models"
HIGH_N = 10            # "high-n_i" threshold for the corr(d,u) cut
CI = (2.5, 97.5)
POPS = ("ALL", "Po10")


# --------------------------------------------------------------------------- #
# locating the production AxD (full) fit                                       #
# --------------------------------------------------------------------------- #
def full_dir(spec, nu: float) -> Path:
    cfg = drift_cfg(nu, variant="full")
    parent = MODELS_ROOT / registry.slice_slug(spec)
    return registry.fit_path(parent, drift_stem("full", nu), spec, cfg, resample_tag="base")


def present(fit_dir: Path) -> bool:
    return (fit_dir / "fit.pkl").is_file() and (fit_dir / "manifest.json").is_file()


def _an_tilde(fd) -> np.ndarray:
    """Within-athlete-centered career age (length N) = A_n - mean_i(A_n)."""
    row = fd.row_idx
    n = np.bincount(row, minlength=fd.I)
    s = np.bincount(row, weights=fd.A_n, minlength=fd.I)
    bar = np.where(n > 0, s / np.maximum(n, 1.0), 0.0)
    return fd.A_n - bar[row]


# --------------------------------------------------------------------------- #
# the loaded fit + rebuilt EB posterior var                                   #
# --------------------------------------------------------------------------- #
@dataclass
class DriftFit:
    cohort: str
    sex: str
    slug: str
    nutag: str
    fd: object
    elig: np.ndarray      # (I,) bool
    d: np.ndarray         # (I,) EB drift slope (0 where ineligible)
    u: np.ndarray         # (I,) athlete factor
    n_i: np.ndarray       # (I,) race count
    A_e: np.ndarray       # (I,) entry age (age at debut)
    post_var: np.ndarray  # (I,) EB posterior var (nan where ineligible)
    omega_d2: float
    model: object = None  # live fitted model (for aging-curve reconstruction)

    def table(self, min_n: int = 3) -> pd.DataFrame:
        """Eligible athletes with n_i >= min_n."""
        e = self.elig & (self.n_i >= min_n)
        return pd.DataFrame({
            "u": self.u[e], "d": self.d[e], "n_i": self.n_i[e].astype(int),
            "A_e": self.A_e[e], "post_var": self.post_var[e],
            "improver": self.d[e] < 0,
        }).reset_index(drop=True)


def load_drift(cohort: str, sex: str, *, nu: float = 8.0, mrc: int = 2) -> DriftFit | None:
    name = f"{cohort}_{sex}"
    spec = S.build_spec(name, min_race_count=mrc)
    slug = registry.slice_slug(spec)
    fdir = full_dir(spec, nu)
    if not present(fdir):
        print(f"  [skip] {slug}: no full fit at {fdir}")
        return None

    fd = load_slice(spec)
    m = registry.load_fit(fdir, fd)
    elig = np.asarray(m.eligible_d, dtype=bool)
    d = np.asarray(m.params["d"], dtype=float)
    u = np.asarray(m.params["u"], dtype=float)
    omega_d2 = float(m.params.get("omega_d2", np.nan))
    n_i = np.bincount(fd.row_idx, minlength=fd.I).astype(float)
    A_e = np.asarray(fd.A_e, dtype=float)

    # rebuild the d-block Fisher info s_den_i = sum w_base * ac^2  (post_var unsaved)
    sigma2 = float(m.params.get("sigma2", np.nan))
    nu_val = float(m.params.get("nu", nu))
    r = m.residuals()
    w_base = irls_weights(r, sigma2, nu_val) / max(sigma2, 1e-300)
    ac = _an_tilde(fd)
    s_den = np.bincount(fd.row_idx, weights=w_base * ac * ac, minlength=fd.I)
    ridge = 1.0 / omega_d2 if omega_d2 > 0 else 0.0
    post_var = np.where(elig, 1.0 / np.maximum(s_den + ridge, 1e-300), np.nan)

    nutag = f"nu{nu:.2f}".replace(".", "p")
    return DriftFit(cohort=cohort, sex=sex, slug=slug, nutag=nutag, fd=fd,
                    elig=elig, d=d, u=u, n_i=n_i, A_e=A_e, post_var=post_var,
                    omega_d2=omega_d2, model=m)


@functools.lru_cache(maxsize=None)
def load_cached(cohort: str, sex: str, nu: float, mrc: int) -> DriftFit | None:
    return load_drift(cohort, sex, nu=nu, mrc=mrc)


# --------------------------------------------------------------------------- #
# descriptors                                                                 #
# --------------------------------------------------------------------------- #
def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    """Fast Spearman = Pearson of ranks (no p-value; called many times in boot)."""
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return np.nan
    rx = x[m].argsort().argsort()
    ry = y[m].argsort().argsort()
    if rx.std() < 1e-12 or ry.std() < 1e-12:
        return np.nan
    return float(np.corrcoef(rx, ry)[0, 1])


def disatt_var(d: np.ndarray, post_var: np.ndarray) -> float:
    """True drift variance = Var(d_hat) - mean(EB posterior var). Clamped >= 0."""
    return max(float(np.var(d, ddof=1) - np.nanmean(post_var)), 0.0)


def descriptors(d: np.ndarray, u: np.ndarray, post_var: np.ndarray,
                n_i: np.ndarray, *, high_n: int = HIGH_N) -> dict:
    hi = n_i >= high_n
    return dict(
        n=int(len(d)),
        disatt_sd=float(np.sqrt(disatt_var(d, post_var))),
        raw_sd=float(np.std(d, ddof=1)),
        skew=float(stats.skew(d)),
        frac_improver=float(np.mean(d < 0)),
        corr_du=_spearman(d, u),
        corr_du_highn=(_spearman(d[hi], u[hi]) if hi.sum() >= 3 else np.nan),
        n_highn=int(hi.sum()),
    )


# descriptor -> (label, gauge-invariant?)
DESC_META = {
    "disatt_sd":     ("disatt SD (log/yr)", True),
    "skew":          ("skew", True),
    "frac_improver": ("frac improver", False),
    "corr_du":       ("corr(d,u)", False),
    "corr_du_highn": ("corr(d,u) hi-n", False),
}
BOOT_KEYS = list(DESC_META)


def bootstrap(d, u, post_var, n_i, *, B: int = 2000, seed: int = 0) -> dict:
    """Nonparametric athlete bootstrap -> {key: np.ndarray(B)} of resampled draws."""
    Ne = len(d)
    rng = np.random.default_rng(seed)
    out = {k: np.empty(B) for k in BOOT_KEYS}
    for b in range(B):
        idx = rng.integers(0, Ne, Ne)
        dd = descriptors(d[idx], u[idx], post_var[idx], n_i[idx])
        for k in BOOT_KEYS:
            out[k][b] = dd[k]
    return out
