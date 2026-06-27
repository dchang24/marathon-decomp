"""Shape metrics for a reconstructed aging curve (pure numpy, no fitting).

All metrics are defined on the **debut-centered** curve ``c(A_n) = g(A_n) - g(0)``
where ``g(A_n) = theta_aging @ B(A_n)`` (optionally plus the entry-age gamma block
at a chosen entry age). Centering on the debut value makes ``A_n = 0`` the zero
reference, so "performance better than first race" is exactly the region
``c(A_n) < 0`` (log-time is lower => faster).

The grid ``A_grid`` is assumed to start at ``A_n = 0`` and be uniformly spaced
(as built by ``e01_aging_grid``). Lengths are reported in years of career age.

These are recomputed post-fit by ``q01_curve_metrics`` from the persisted fit
payloads -- they never trigger a re-fit.
"""
from __future__ import annotations

import numpy as np


def _first_upcross(A: np.ndarray, c: np.ndarray, start: int) -> float:
    """First A where c crosses from <0 to >=0 at or after index `start` (linear
    interp between grid points). NaN if it never crosses back up."""
    for k in range(max(start, 1), len(c)):
        if c[k - 1] < 0.0 <= c[k]:
            # linear interpolation of the zero crossing
            x0, x1, y0, y1 = A[k - 1], A[k], c[k - 1], c[k]
            if y1 == y0:
                return float(x1)
            return float(x0 + (0.0 - y0) * (x1 - x0) / (y1 - y0))
    return float("nan")


def curve_shape_metrics(
    g: np.ndarray,
    A_grid: np.ndarray,
    *,
    entry_age: float | None = None,
    an_p95: float | None = None,
) -> dict:
    """Flat dict of shape metrics for one curve ``g`` on ``A_grid``.

    Parameters
    ----------
    g : the raw aging curve values (log-time units) on A_grid.
    entry_age : athlete entry age (years); when given, ``peak_age`` is reported
        as ``entry_age + peak_an`` so the peak is on a calendar-age axis.
    an_p95 : the 95th percentile of *observed* A_n in the slice; the tail
        wiggle metrics are computed on ``A_n >= an_p95`` (extrapolation zone).
        Falls back to the 95th percentile of the grid when None.
    """
    A = np.asarray(A_grid, dtype=np.float64)
    g = np.asarray(g, dtype=np.float64)
    c = g - g[0]                                   # debut-centered
    dx = float(A[1] - A[0]) if len(A) > 1 else 0.0

    peak_idx = int(np.argmin(c))
    peak_an = float(A[peak_idx])
    peak_depth = float(c[peak_idx])                # <= 0
    peak_pct = float(np.expm1(peak_depth))         # fractional speed-up (negative)

    out: dict[str, float] = {
        "peak_an": peak_an,
        "peak_age": float(entry_age + peak_an) if entry_age is not None
        and np.isfinite(entry_age) else float("nan"),
        "peak_depth": peak_depth,
        "peak_pct": peak_pct,
        "improve_span_yr": float(dx * np.sum(c < 0.0)),
        "breakeven_an": _first_upcross(A, c, peak_idx),
        "tail_value": float(c[-1]),
    }

    # plateau widths: A_n length where the curve is within p% of the peak depth.
    for p in (0.05, 0.10, 0.25):
        thresh = (1.0 - p) * peak_depth            # less-negative bound
        out[f"plateau_p{int(p * 100):02d}_yr"] = float(dx * np.sum(c <= thresh))

    # tail wiggle: sign changes of the slope beyond the observed-A_n p95.
    cut = float(an_p95) if an_p95 is not None and np.isfinite(an_p95) \
        else float(np.percentile(A, 95))
    tail = A >= cut
    if tail.sum() >= 3:
        dc = np.diff(c[tail])
        sign = np.sign(dc[dc != 0.0])
        sign_changes = int(np.sum(sign[1:] != sign[:-1])) if sign.size > 1 else 0
    else:
        sign_changes = 0
    out["tail_an_cut"] = cut
    out["tail_sign_changes"] = float(sign_changes)
    out["tail_nonmonotone"] = float(sign_changes > 0)
    return out
