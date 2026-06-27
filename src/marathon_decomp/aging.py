"""Aging-basis construction and curve reconstruction.

Model-agnostic pure functions used by the model's aging block:

  basis_kind : {"poly", "spline"}
      "poly":   {A_n, A_n^2, ..., A_n^degree}, optionally QR-orthogonalized.
      "spline": natural cubic regression spline (ESL II 5.2.1).

  gamma_form : {"scalar", "varying"}
      "scalar":  entry-age block = gamma * (A_e - mean) * A_n           (1 col)
      "varying": entry-age block = sum_k gamma_k * (A_e - mean) * B_k   (K cols)

The fitted aging coefficients are `theta_aging` (the phi_k, one per basis
column) and `gamma` (the gamma_k, length 1 for scalar / K for varying).

Curve reconstruction comes in two flavours:
  - `*_on_grid(model, ...)`  — from a live fitted model.
  - `*_from_payload(payload, ...)` — from a saved pickle dict alone, using the
    basis metadata persisted in `payload["model_extra"]`. The payload form is
    self-contained: it needs neither the model object nor the original FitData.
"""
from __future__ import annotations

from typing import Any

import numpy as np


# ---------------------------------------------------------------------------
# Basis builders (pure functions of A_n)
# ---------------------------------------------------------------------------

def _raw_poly_columns(A_n: np.ndarray, degree: int) -> np.ndarray:
    """Raw monomial design matrix {A_n, A_n^2, ..., A_n^degree}, shape (N, degree).

    No constant column — the global intercept is absorbed into u_i.
    """
    if degree < 1:
        raise ValueError(f"poly degree must be >= 1, got {degree}")
    cols = [A_n.astype(np.float64).copy()]
    for k in range(2, degree + 1):
        cols.append(A_n.astype(np.float64) ** k)
    return np.column_stack(cols)


def build_poly_basis(
    A_n: np.ndarray, degree: int, orthogonalize: bool = True,
) -> tuple[np.ndarray, tuple]:
    """Polynomial basis on A_n. Returns (B, transform).

    `transform` is a (R, scale) pair that lets the fitted theta be applied
    to any other A_n grid via `poly_curve_from_theta`. When
    `orthogonalize=True`, columns of B are orthonormal then rescaled by
    sqrt(N) so all entries are O(1) — this kills the inner-block
    correlation that makes raw monomials slow under CD.
    """
    B_raw = _raw_poly_columns(A_n, degree)
    if not orthogonalize:
        return B_raw, (np.eye(degree), 1.0)
    Q, R = np.linalg.qr(B_raw)
    scale = float(np.sqrt(B_raw.shape[0]))
    return Q * scale, (R, scale)


def poly_curve_from_theta(
    theta: np.ndarray, transform: tuple, A_grid: np.ndarray, degree: int,
) -> np.ndarray:
    """Reconstruct the fitted aging curve on an arbitrary A_n grid.

    Inverts the QR + rescale used at fit time.
    """
    R, scale = transform
    B_raw = _raw_poly_columns(A_grid, degree)
    if scale == 1.0 and np.allclose(R, np.eye(degree)):
        B = B_raw
    else:
        # Q_grid = B_raw @ R^{-1}, then rescale.
        B = np.linalg.solve(R.T, B_raw.T).T * scale
    return B @ theta


def ncs_basis(x: np.ndarray, knots: np.ndarray) -> np.ndarray:
    """Natural cubic regression spline basis on x with given knots.

    Returns (N, K-1) matrix: 1 linear column + K-2 cubic-correction columns
    (no constant column — absorbed by u_i). See ESL II 5.2.1.
    """
    knots = np.asarray(knots, dtype=np.float64)
    K = len(knots)
    if K < 3:
        raise ValueError(f"need >= 3 knots, got {K}")
    xK = knots[-1]
    xv = x.astype(np.float64)

    def d_fn(k: int) -> np.ndarray:
        num1 = np.maximum(xv - knots[k], 0.0) ** 3
        num2 = np.maximum(xv - xK, 0.0) ** 3
        return (num1 - num2) / (xK - knots[k])

    cols = [xv.copy()]
    d_last = d_fn(K - 2)
    for k in range(K - 2):
        cols.append(d_fn(k) - d_last)
    return np.column_stack(cols)


def default_knots_from_An(
    A_n: np.ndarray, n_knots: int = 5, q_lo: float = 0.05, q_hi: float = 0.95,
) -> np.ndarray:
    """Default spline knots: 0 (debut anchor) + interior quantiles of positive A_n.

    A_n has a point mass at 0 (every athlete's debut race) and a sparse upper
    tail (a handful of athletes racing >10 yr after debut). Two choices keep
    the knots well-placed:

      * Quantiles are taken on the **positive subset only**, so the debut mass
        does not drag the interior knots toward 0.
      * The ``n_knots - 1`` positive knots span the **[q_lo, q_hi] quantile
        range** (default p5..p95), not [min, max]. This avoids both (a) a
        near-duplicate knot at the smallest positive A_n (~0.003 yr, the old
        0.0-quantile) and (b) a boundary knot pinned to the lone maximum out in
        the sparse tail (the old 1.0-quantile). The explicit 0 is kept as the
        lower boundary (A_n cannot be negative).
    """
    A = np.asarray(A_n, dtype=np.float64)
    A = A[np.isfinite(A)]
    pos = A[A > 0]
    if len(pos) == 0:
        return np.linspace(0.0, 1.0, n_knots)
    qs = np.linspace(q_lo, q_hi, n_knots - 1)
    inner = np.quantile(pos, qs)
    knots = np.concatenate([[0.0], inner])
    # ensure strictly increasing
    for i in range(1, len(knots)):
        if knots[i] <= knots[i - 1]:
            knots[i] = knots[i - 1] + 1e-3
    return knots


# ---------------------------------------------------------------------------
# Curve reconstruction — from a live model
# ---------------------------------------------------------------------------

def aging_curve_on_grid(model: Any, A_grid: np.ndarray) -> np.ndarray:
    """Reconstruct the fitted aging curve theta_aging @ B(A_grid).

    For poly bases this applies the same QR transform used at fit time; for
    spline bases it evaluates the natural cubic spline on the stored knots.
    The returned curve does NOT include the entry-age gamma term (which is
    athlete-specific via A_e); see `entry_age_curve_on_grid`.
    """
    cfg = model.config
    theta = model.params["theta_aging"]
    if cfg.basis_kind == "poly":
        return poly_curve_from_theta(theta, model._basis_transform, A_grid, cfg.degree)
    if cfg.basis_kind == "spline":
        return ncs_basis(A_grid, model._spline_knots) @ theta
    raise ValueError(f"unknown basis_kind={cfg.basis_kind!r}")


def entry_age_curve_on_grid(
    model: Any, A_grid: np.ndarray, A_e_centered: float
    ) -> np.ndarray:
    """Reconstruct gamma_block(A_grid) for a given centred entry age.

    For scalar gamma:   gamma * A_e_centered * A_grid.
    For varying gamma:  A_e_centered * sum_k gamma_k * B_k(A_grid).
    """
    cfg = model.config
    gamma = np.atleast_1d(model.params["gamma"])
    if cfg.gamma_form == "scalar":
        return float(gamma[0]) * A_e_centered * np.asarray(A_grid, dtype=np.float64)
    # varying gamma — same basis on A_grid as the aging block.
    if cfg.basis_kind == "poly":
        B = _poly_design_on_grid(A_grid, model._basis_transform, cfg.degree)
    else:
        B = ncs_basis(np.asarray(A_grid, dtype=np.float64), model._spline_knots)
    return A_e_centered * (B @ gamma)


# ---------------------------------------------------------------------------
# Curve reconstruction — from a saved payload (self-contained)
# ---------------------------------------------------------------------------

def aging_curve_from_payload(payload: dict, A_grid: np.ndarray) -> np.ndarray:
    """Reconstruct the aging curve from a `BaseModel.save()` pickle dict alone.

    Uses `payload["model_extra"]` (basis metadata) + `payload["params"]`;
    needs neither the model object nor the original FitData.
    """
    meta = payload["model_extra"]
    theta = np.asarray(payload["params"]["theta_aging"], dtype=np.float64)
    if meta["basis_kind"] == "poly":
        return poly_curve_from_theta(theta, meta["poly_transform"], A_grid, meta["degree"])
    if meta["basis_kind"] == "spline":
        return ncs_basis(A_grid, np.asarray(meta["spline_knots"])) @ theta
    raise ValueError(f"unknown basis_kind={meta['basis_kind']!r}")


def entry_age_curve_from_payload(
    payload: dict, A_grid: np.ndarray, A_e_centered: float,
) -> np.ndarray:
    """Reconstruct gamma_block(A_grid) from a saved pickle dict alone."""
    meta = payload["model_extra"]
    gamma = np.atleast_1d(np.asarray(payload["params"]["gamma"], dtype=np.float64))
    if meta["gamma_form"] == "scalar":
        return float(gamma[0]) * A_e_centered * np.asarray(A_grid, dtype=np.float64)
    if meta["basis_kind"] == "poly":
        B = _poly_design_on_grid(A_grid, meta["poly_transform"], meta["degree"])
    else:
        B = ncs_basis(np.asarray(A_grid, dtype=np.float64), np.asarray(meta["spline_knots"]))
    return A_e_centered * (B @ gamma)


def _poly_design_on_grid(A_grid: np.ndarray, transform: tuple, degree: int) -> np.ndarray:
    """Poly design matrix on A_grid, matching the fit-time QR reparametrization."""
    R, scale = transform
    B_raw = _raw_poly_columns(np.asarray(A_grid, dtype=np.float64), degree)
    if scale == 1.0 and np.allclose(R, np.eye(degree)):
        return B_raw
    return np.linalg.solve(R.T, B_raw.T).T * scale
