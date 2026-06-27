"""Resampling weights for uncertainty quantification of the latent factors.

Currently implements the **Bayesian (weighted) bootstrap over athletes** — the
recommended scheme for quantifying how the race-side and global factors wobble
as the athlete sample varies.

Why Bayesian, not the classic with-replacement athlete bootstrap: the design is
gatekept by `min_runners_per_race` and reduced to a single connected component
(see `data.load_slice`). A nonparametric bootstrap omits ~37% of athletes per
replicate, which drops races below the field-size floor and can fragment the
graph — so each replicate has a *different* race support and `v_j` is missing
for many races. Strictly-positive Dirichlet weights keep every athlete (and
hence every race) in every replicate, so the support is fixed at the full-data
design and per-race intervals are clean. The weights are also a smooth
perturbation of unit weights, so warm-starting each replicate from the
full-data fit converges in a few iterations.

How it plugs in: the model is weighted least squares. A per-athlete weight
`g_i` expanded to cells (`boot_cell_weights`) is set on the fitter via
`model.boot_w_cell` and folded into the IRLS weight, so it flows uniformly into
every weighted block (u, v, d, the aging/gamma block, sigma^2 and the EB
variance update). It composes with — does not replace — the runner-reliability
weight `w_ath_i = n_i/(n_i+n0)`: the effective cross-athlete weight becomes
`g_i * w_ath_i * irls/sigma^2`, i.e. (resampling draw) x (reliability).

Note on interpretation: because a single athlete's cells all share one scalar
`g_i`, that scalar cancels in the per-athlete `u_i` update — the Bayesian
*athlete* bootstrap therefore injects almost no direct uncertainty into the
conditional `u_i` estimate (only indirect, via the shifted race factors). It is
the right tool for v / sigma^2 / omega / aging uncertainty, NOT for genuine
per-athlete ability sampling error (that needs a within-athlete cell bootstrap).

Public API:
    bayesian_athlete_weights(fd, rng, concentration=1.0) -> g (length I, mean 1)
    boot_cell_weights(fd, g)                              -> g[row] (length N)
"""
from __future__ import annotations

import numpy as np

from .data import FitData

__all__ = ["bayesian_athlete_weights", "boot_cell_weights"]


def bayesian_athlete_weights(
    fd: FitData,
    rng: np.random.Generator,
    *,
    concentration: float = 1.0,
) -> np.ndarray:
    """Dirichlet bootstrap weights over the I athletes, scaled to mean 1.

    Draws `g ~ Dirichlet(concentration, ..., I times)` (via Gamma draws) and
    rescales so the sample mean is exactly 1, hence `sum_i g_i = I`. This is the
    smooth analog of a multinomial athlete bootstrap whose counts sum to I: in
    expectation the total cell weight `sum_i g_i * n_i = N`, so the model's
    `sigma2 = sum(w * r^2) / N` stays correctly normalized.

    `concentration=1.0` is the standard Rubin Bayesian bootstrap (flat
    Dirichlet). Larger values shrink the weights toward uniform (less
    resampling variance); smaller values inflate it.

    Returns a float64 array of length `fd.I`.
    """
    if concentration <= 0:
        raise ValueError(f"concentration must be > 0, got {concentration}")
    g = rng.gamma(shape=concentration, scale=1.0, size=fd.I)
    m = g.mean()
    if m <= 0:  # astronomically unlikely; guard anyway
        return np.ones(fd.I, dtype=np.float64)
    return (g / m).astype(np.float64)


def boot_cell_weights(fd: FitData, g: np.ndarray) -> np.ndarray:
    """Expand per-athlete weights `g` (length I) to per-cell weights (length N)
    by indexing with `fd.row_idx`. The result is what the fitter consumes via
    `model.boot_w_cell`."""
    g = np.asarray(g, dtype=np.float64)
    if g.shape != (fd.I,):
        raise ValueError(f"g must have shape ({fd.I},), got {g.shape}")
    return g[fd.row_idx]
