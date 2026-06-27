"""K-fold cross-validation primitives for the latent-factor models.

The model only identifies an athlete factor u_i if athlete i has training data,
so naive per-observation folds are wrong: an athlete with all finishes in the
held-out fold becomes unpredictable. We stratify **within athlete** — deal each
athlete's finishes round-robin into K folds — so every athlete keeps at least
``n_i - ceil(n_i/K) >= 1`` training finishes in every fold (no singletons), while
each finish is still held out exactly once. Races keep ~(K-1)/K of their field,
so race factors stay identified too.

Public API:
    assign_folds(fd, K, seed, min_test_n)  -> fold id per cell (-1 = always train)
    subset_fitdata(fd, keep_mask)          -> a re-indexed FitData on those cells
    heldout_logdensity(fd, test_mask, m)   -> predictive log-density of held-out cells

Scope: the held-out predictor is u_i + v_j (the rank-1 baseline) **plus the
global aging block** (theta_aging @ B(A_n) and the gamma entry-age term) when
the trained model has them on — the aging regressors are rebuilt on the
held-out cells from the train-fitted basis. The per-athlete drift d_i is NOT
scored out-of-sample (it is an individual term with no held-out value); a model
with use_d on contributes nothing from d_i to the held-out prediction.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.special import gammaln

from .aging import aging_curve_on_grid, entry_age_curve_on_grid
from .base import BaseModel
from .data import FitData


# ---------------------------------------------------------------------------
# Student-t / Gaussian log predictive density
# ---------------------------------------------------------------------------

def t_logpdf(r: np.ndarray, sigma2: float, nu: float) -> np.ndarray:
    """Per-cell log density of residual `r` under t_nu(0, sigma2) (Gaussian if
    nu=inf). Matches BaseModel.log_lik's kernel so CV scores are on the same
    scale as the in-sample loglik."""
    s2 = max(float(sigma2), 1e-300)
    if not np.isfinite(nu):
        return -0.5 * (np.log(2.0 * np.pi * s2) + r * r / s2)
    return (
        gammaln((nu + 1) / 2) - gammaln(nu / 2) - 0.5 * np.log(nu * np.pi * s2)
        - 0.5 * (nu + 1) * np.log1p(r * r / (nu * s2))
    )


# ---------------------------------------------------------------------------
# Fold assignment
# ---------------------------------------------------------------------------

def assign_folds(
    fd: FitData, K: int, seed: int = 0, min_test_n: int = 1,
) -> np.ndarray:
    """Stratified-within-athlete fold id per cell (length N), in 0..K-1.

    Each athlete's cells are shuffled and dealt round-robin into the K folds
    (with a random per-athlete rotation so small-n athletes don't all seed
    fold 0). Cells of athletes with fewer than `min_test_n` finishes are set to
    -1 ("always train, never tested") — `min_test_n=1` (default) tests everyone.
    """
    if K < 2:
        raise ValueError(f"K must be >= 2, got {K}")
    rng = np.random.default_rng(seed)
    N, I = fd.N, fd.I
    row = fd.row_idx

    # sort cells by athlete, random within athlete
    order = np.lexsort((rng.random(N), row))
    row_sorted = row[order]
    _, first_idx, counts = np.unique(row_sorted, return_index=True, return_counts=True)
    rank = np.arange(N) - np.repeat(first_idx, counts)        # position within athlete
    off = rng.integers(0, K, size=I)                          # per-athlete rotation
    fold_sorted = (rank + off[row_sorted]) % K
    fold = np.empty(N, dtype=np.int64)
    fold[order] = fold_sorted

    if min_test_n > 1:
        n_i = np.bincount(row, minlength=I)
        fold[n_i[row] < min_test_n] = -1
    return fold


# ---------------------------------------------------------------------------
# Build a re-indexed FitData on a subset of cells
# ---------------------------------------------------------------------------

def subset_fitdata(fd: FitData, keep_mask: np.ndarray) -> FitData:
    """A new FitData containing only the kept cells, reduced to its largest
    connected component and re-indexed 0..I'-1 / 0..J'-1.

    `athlete_ids` / `race_ids` keep the ORIGINAL surrogate ids, so a fit on the
    subset can be matched back to held-out cells by id. Per-cell A_n and
    per-athlete A_e are carried from the parent (correct for the baseline, which
    ignores them; recompute for aging models).
    """
    orig = np.flatnonzero(np.asarray(keep_mask, dtype=bool))
    if orig.size == 0:
        raise ValueError("keep_mask selects no cells")
    row, col = fd.row_idx[orig], fd.col_idx[orig]

    # restrict to the giant component of the kept bipartite graph
    _, ri = np.unique(row, return_inverse=True)
    uc, ci = np.unique(col, return_inverse=True)
    A, R = ri.max() + 1, len(uc)
    g = coo_matrix((np.ones(len(ri)), (ri, A + ci)), shape=(A + R, A + R))
    n_comp, labels = connected_components(g, directed=False)
    if n_comp > 1:
        cell_label = labels[ri]
        giant = np.bincount(cell_label).argmax()
        gk = cell_label == giant
        orig, row, col = orig[gk], row[gk], col[gk]

    # original ids per kept cell, then re-factorize to dense 0..I'-1 / 0..J'-1
    aid = fd.athlete_ids[row]
    rid = fd.race_ids[col]
    arow, auniq = pd.factorize(aid, sort=True)
    acol, runiq = pd.factorize(rid, sort=True)

    # map surviving ids back to parent per-athlete / per-race rows
    aidx = {int(a): i for i, a in enumerate(fd.athlete_ids)}
    ridx = {int(r): j for j, r in enumerate(fd.race_ids)}
    ai = np.array([aidx[int(a)] for a in auniq], dtype=np.int64)
    rj = np.array([ridx[int(r)] for r in runiq], dtype=np.int64)

    return FitData(
        row_idx=arow.astype(np.int64),
        col_idx=acol.astype(np.int64),
        y=fd.y[orig].copy(),
        A_n=fd.A_n[orig].copy(),
        A_e=fd.A_e[ai].copy(),
        athlete_ids=fd.athlete_ids[ai].copy(),
        athlete_sex=fd.athlete_sex[ai].copy(),
        athlete_country=fd.athlete_country[ai].copy(),
        race_ids=fd.race_ids[rj].copy(),
        race_series=fd.race_series[rj].copy(),
        race_country=fd.race_country[rj].copy(),
        race_date=fd.race_date[rj].copy(),
        I=len(auniq), J=len(runiq), N=len(orig),
        spec=fd.spec,
        data_version=fd.data_version,
        cache_key=f"{fd.cache_key}_cv{len(orig)}",
        response_kind=fd.response_kind,
    )


# ---------------------------------------------------------------------------
# Held-out scoring
# ---------------------------------------------------------------------------

def heldout_logdensity(
    fd: FitData, test_mask: np.ndarray, train_model: BaseModel,
) -> dict:
    """Predictive log-density of the held-out cells under a train-fitted model.

    For each test cell (i, j) the prediction is u_i + v_j looked up from the
    train fit by ORIGINAL athlete/race id (so the train re-indexing is handled),
    **plus the train-fitted aging block** theta@B(A_n) and the gamma entry-age
    term when the model has them on (the basis is re-evaluated at the held-out
    cells' A_n; A_e is centered by the train mean tm.Ae_bar, NaN-A_e -> 0 as at
    fit time). Scored against the train-fitted (sigma2, nu). Test cells whose
    athlete or race is absent from the train fit (orphans dropped by the
    giant-component reduction) are excluded and counted.

    Returns: sum_logdens, mean_logdens, n_test, n_orphan, rmse.
    """
    test = np.flatnonzero(np.asarray(test_mask, dtype=bool))
    aid = pd.Series(fd.athlete_ids[fd.row_idx[test]])
    rid = pd.Series(fd.race_ids[fd.col_idx[test]])
    y = fd.y[test]

    tm = train_model
    u_by_aid = pd.Series(tm.params["u"], index=tm.data.athlete_ids)
    v_by_rid = pd.Series(tm.params["v"], index=tm.data.race_ids)
    u = aid.map(u_by_aid).to_numpy()
    v = rid.map(v_by_rid).to_numpy()

    pred = u + v
    cfg = tm.config
    if getattr(cfg, "use_phi12", False) and getattr(tm, "K_basis", 0) > 0:
        # theta @ B(A_n) re-evaluated on the held-out cells via the train basis.
        pred = pred + aging_curve_on_grid(tm, fd.A_n[test])
    if getattr(cfg, "use_gamma", False):
        ae_test = fd.A_e[fd.row_idx[test]]
        ae_c = np.where(np.isnan(ae_test), 0.0, ae_test - float(getattr(tm, "Ae_bar", 0.0)))
        pred = pred + entry_age_curve_on_grid(tm, fd.A_n[test], ae_c)

    ok = np.isfinite(u) & np.isfinite(v)
    n_orphan = int((~ok).sum())
    r = y[ok] - pred[ok]

    sigma2 = float(tm.params["sigma2"])
    nu = float(tm.params["nu"])
    logd = t_logpdf(r, sigma2, nu)
    return {
        "sum_logdens": float(logd.sum()),
        "mean_logdens": float(logd.mean()) if r.size else float("nan"),
        "n_test": int(r.size),
        "n_orphan": n_orphan,
        "rmse": float(np.sqrt(np.mean(r * r))) if r.size else float("nan"),
    }
