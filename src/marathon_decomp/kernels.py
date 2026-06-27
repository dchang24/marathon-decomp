"""Numba-compiled and small numpy primitives shared by all model fitters.

No model semantics live here. Each function is a numerical primitive:
inner-CD sweep over a few correlated scalar regressors, IRLS reweighting,
etc. Model files in `models/` compose these into a full outer loop.
"""
from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True, fastmath=True)
def cd_inner_scalars(
    r: np.ndarray,
    w: np.ndarray,
    x_cols: np.ndarray,
    use_mask: np.ndarray,
    theta: np.ndarray,
    n_iters: int,
) -> None:
    """In-place block coordinate descent over `theta`'s active entries.

    Solves the small correlated-regressor sub-problem (the aging/gamma
    coefficients) by sweeping one variable at a time, patching the running
    residual `r` between scalar updates so each next coordinate sees a
    fresh residual. After `n_iters` full sweeps the active entries of
    `theta` are at the joint conditional minimum w.r.t. the supplied
    weights, and `r` is consistent with the final values.

    Args:
        r: length-N running residual; modified in place.
        w: length-N per-cell weight (callers form this — typically
            irls_w/sigma^2 * w_ath).
        x_cols: shape (K, N) regressor matrix.
        use_mask: shape (K,) bool — skip k where False.
        theta: shape (K,) current values; modified in place.
        n_iters: number of full sweeps over the K coordinates.
    """
    K = theta.shape[0]
    N = r.shape[0]
    eps = 1e-10
    for _ in range(n_iters):
        for k in range(K):
            if not use_mask[k]:
                continue
            old = theta[k]
            num = 0.0
            den = 0.0
            for i in range(N):
                xi = x_cols[k, i]
                wxi = w[i] * xi
                # running-residual identity: partial_resid = r + old*x
                num += wxi * (r[i] + old * xi)
                den += wxi * xi
            new = num / (den + eps)
            delta = new - old
            theta[k] = new
            if delta != 0.0:
                for i in range(N):
                    r[i] -= delta * x_cols[k, i]


def irls_weights(r: np.ndarray, sigma2: float, nu: float) -> np.ndarray:
    """Per-cell Student-t IRLS weights w_ij = (nu+1) / (nu + r^2/sigma^2).

    Returns 1.0 everywhere when `nu = inf` so the Gaussian limit recovers
    plain L2 bit-for-bit (not an evaluation of the t-density at infinity).
    """
    if not np.isfinite(nu):
        return np.ones_like(r)
    return (nu + 1.0) / (nu + (r * r) / max(sigma2, 1e-300))
