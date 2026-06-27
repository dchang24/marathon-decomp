"""Reproducible (u, v) initialization schemes for the rank-1 baseline.

For the baseline only ``u`` and ``v`` are free (s/d/aging off), so a warmstart
dict need only set those two; the fitter zero-fills the rest. Every random init
is fully determined by ``(FitData, seed, jit_u, jit_v)`` so it can be recreated
later from its recorded seed.
"""
from __future__ import annotations

import numpy as np

from marathon_decomp import FitData


def mean_u(fd: FitData) -> np.ndarray:
    """Per-athlete mean response — the canonical 'mean' init for u."""
    sum_y = np.bincount(fd.row_idx, weights=fd.y, minlength=fd.I)
    n_i = np.bincount(fd.row_idx, minlength=fd.I)
    return np.where(n_i > 0, sum_y / np.maximum(n_i, 1.0), 0.0)


def perturb(base: dict[str, np.ndarray], seed: int, *,
            jit_u: float, jit_v: float) -> dict[str, np.ndarray]:
    """Random Gaussian perturbation around an arbitrary (u, v) base.

    Deterministic in `seed`: same base + same seed reproduces it exactly.
    """
    rng = np.random.default_rng(seed)
    u = np.asarray(base["u"], dtype=np.float64)
    v = np.asarray(base["v"], dtype=np.float64)
    return {"u": u + rng.normal(0.0, jit_u, size=u.shape),
            "v": v + rng.normal(0.0, jit_v, size=v.shape)}


def random_init(fd: FitData, seed: int, *,
                jit_u: float, jit_v: float) -> dict[str, np.ndarray]:
    """Random Gaussian perturbation around the mean init."""
    return perturb({"u": mean_u(fd), "v": np.zeros(fd.J)},
                   seed, jit_u=jit_u, jit_v=jit_v)


def build_inits(fd: FitData, *, n_random: int, jit_u: float, jit_v: float,
                seed0: int, base: dict[str, np.ndarray] | None = None,
                ) -> list[tuple[str, int | None, dict[str, np.ndarray]]]:
    """Return ``[(name, seed, warmstart), ...]``.

    A deterministic anchor init (seed=None) followed by `n_random` perturbed
    restarts with seeds ``seed0 .. seed0+n_random-1``. When `base` is given the
    anchor is that (u, v) and restarts perturb around it (warm continuation,
    e.g. the L2 solution); otherwise the anchor is the mean init named ``mean``.
    """
    if base is None:
        anchor_name = "mean"
        anchor = {"u": mean_u(fd), "v": np.zeros(fd.J)}
    else:
        anchor_name = "warm"
        anchor = {"u": np.asarray(base["u"]).copy(), "v": np.asarray(base["v"]).copy()}
    inits: list[tuple[str, int | None, dict[str, np.ndarray]]] = [
        (anchor_name, None, anchor),
    ]
    for k in range(n_random):
        seed = seed0 + k
        inits.append((f"rand{k}", seed,
                      perturb(anchor, seed, jit_u=jit_u, jit_v=jit_v)))
    return inits
