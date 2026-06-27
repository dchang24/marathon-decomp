"""Shared loaders + rank machinery for the race_comparison scripts.

Library module (no eNN/qNN/pNN prefix): not runnable on its own.

Provides:
  find_fit_dir      locate `{model}_{nutag}_best__{hash}` under results/models/{slug}
  load_editions     point fit -> per-race DataFrame (race_id, v, n_j, series_key,
                    country, date, year) + the fit dir
  load_boot_wide    bootstrap race factors -> (R, n_races) array aligned to a
                    race_id order (replicates share the fit's gauge, see
                    e10_bootstrap docstring -- directly comparable across runs)
  rank_matrix       per-replicate ranks (1 = largest of the row)
  pct               percentile helper

Conventions used by q02/q03:
  "slowest" = largest v_j (rank 1 = slowest), "fastest" = smallest v_j.
  All v are in log-time, beta=0 gauge; x180 ~= minutes at a 3:00 marathon.
"""
from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from marathon_decomp import load_slice  # noqa: E402
from marathon_decomp.config import DATA_DIR, RESULTS_DIR  # noqa: E402

MODELS_ROOT = RESULTS_DIR / "models"
OUT_ROOT = RESULTS_DIR / "analysis" / "race_comparison"
REF_MARATHON_MIN = 180.0   # minutes at a 3:00 marathon


def v_to_min(v):
    """Convert a log-time race factor v (a *level*) to minutes for a 3:00 runner.

    v is a log-time effect, so the exact minute difference vs the average race is
    ``180 * (exp(v) - 1)``, NOT the first-order ``180 * v``: the two diverge for
    large |v| (~1 min at v=0.11). Use this for per-race/per-series difficulty
    levels. It is *not* appropriate for log *contrasts* (a difference of two v's,
    e.g. cross-slice movers or M-vs-W), where x180 stays the natural linearization.
    """
    return REF_MARATHON_MIN * (np.exp(v) - 1.0)


def find_fit_dir(slug: str, model: str, nutag: str) -> Path | None:
    """First `{model}_{nutag}_best__{hash}` dir under results/models/{slug} with a fit.pkl."""
    sd = MODELS_ROOT / slug
    if not sd.is_dir():
        return None
    cands = sorted(sd.glob(f"{model}_{nutag}_best__*"))
    cands = [c for c in cands if (c / "fit.pkl").is_file()]
    return cands[0] if cands else None


def load_editions(slug: str, model: str, nutag: str,
                  data_version: str) -> tuple[pd.DataFrame, Path] | None:
    """Point fit -> one row per race: race_id, v, n_j (+ competitions metadata)."""
    fit_dir = find_fit_dir(slug, model, nutag)
    if fit_dir is None:
        print(f"  [skip] no point fit for {slug} ({model}_{nutag})")
        return None
    payload = pickle.load(open(fit_dir / "fit.pkl", "rb"))
    fd = load_slice(payload["spec"], payload.get("data_version", data_version))
    v = np.asarray(payload["params"]["v"], np.float64)
    race_ids = np.asarray(fd.race_ids, np.int64)
    n_j = np.bincount(np.asarray(fd.col_idx, np.int64), minlength=race_ids.size)
    df = pd.DataFrame({"race_id": race_ids, "v_raw": v, "n_j": n_j.astype(int)})
    comps = pd.read_parquet(DATA_DIR / "competitions.parquet",
                            columns=["race_id", "series_key", "country", "date", "year"])
    df = df.merge(comps, on="race_id", how="left")
    df["series_key"] = df["series_key"].fillna("unknown").astype(str)
    df["year"] = df["year"].fillna(0).astype(int)
    # Re-impose the beta=0 APC gauge (fits <= 1.2.0 are off-manifold). `v` is the
    # gauged production estimand; `v_raw` is kept for diagnostics.
    df["v"] = beta0_gauge(df["v_raw"].to_numpy(), frac_years(df["date"].to_numpy()))
    return df, fit_dir


def load_boot_wide(fit_dir: Path, race_ids: np.ndarray,
                   point_v: np.ndarray) -> np.ndarray | None:
    """Bootstrap replicates of v aligned to `race_ids` -> (R, n) or None.

    Each replicate is re-gauged onto its own beta=0 manifold (see `beta0_gauge`),
    consistent with the gauged point v returned by `load_editions` -- replicates
    from fits <= 1.2.0 otherwise carry a wandering APC tilt. `point_v` (used only
    to fill the rare missing race) should be the *gauged* point estimate.

    Races absent from the bootstrap table (should not happen) are filled with
    the point estimate, so they get degenerate CIs without distorting the
    ranks of the others; a warning is printed.
    """
    p = fit_dir / "bootstrap" / "race_factors.parquet"
    if not p.is_file():
        return None
    rf = pd.read_parquet(p, columns=["run_id", "race_id", "v"])
    wide = (rf[rf["run_id"] > 0]
            .pivot(index="run_id", columns="race_id", values="v")
            .reindex(columns=race_ids))
    if wide.empty:
        return None
    B = wide.to_numpy(np.float64)
    bad = ~np.isfinite(B)
    if bad.any():
        print(f"    [warn] {int(bad.any(axis=0).sum())} races missing from "
              f"bootstrap {fit_dir.name}; filled with point v")
        B = np.where(bad, np.broadcast_to(point_v, B.shape), B)
    return beta0_gauge_rows(B, frac_year_for_races(race_ids))


def frac_years(dates) -> np.ndarray:
    """datetime64 array -> fractional years (matches q01 / p01 convention)."""
    rd = pd.DatetimeIndex(pd.to_datetime(dates))
    return (rd.year + (rd.dayofyear - 1) / 365.25).to_numpy(np.float64)


# --------------------------------------------------------------------------- #
# beta=0 APC gauge                                                            #
# --------------------------------------------------------------------------- #
# Fits saved by fitter <= 1.2.0 sit OFF the beta=0 manifold: the fitter's
# per-iteration apply_apc_gauge_beta0() silently no-opped (a _t_j units bug,
# fixed in fitter 1.3.0; see models/model.py), so their v_j carries a
# prior-pinned secular date tilt (slope(v ~ t) != 0). load_editions /
# load_boot_wide therefore re-impose beta=0 at read time so q02/q03/q05/p20
# rank on the production estimand the docstring/paper claim. Idempotent: c ~= 0
# for on-manifold fits (>= 1.3.0), so this is a no-op once fits are refreshed.
def beta0_gauge(v: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Project v onto beta=0: subtract the OLS date slope, re-centre to mean 0.

    This is the production mean_v_zero + beta=0 convention. Prediction-invariant
    (a pure re-gauge), so it is exact on saved fits without refitting.
    """
    v = np.asarray(v, np.float64); t = np.asarray(t, np.float64)
    tc = t - t.mean()
    denom = float(tc @ tc)
    if denom <= 0.0:                       # degenerate (single date) -> level only
        return v - v.mean()
    c = float(((v - v.mean()) @ tc) / denom)
    out = v - c * t
    return out - out.mean()


def beta0_gauge_rows(B: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Row-wise :func:`beta0_gauge` for a (R, n) bootstrap matrix.

    Each replicate is gauged by its OWN slope c_r, so the wandering per-replicate
    APC tilt is removed (it otherwise inflates CIs for early/late-window races).
    """
    B = np.asarray(B, np.float64); t = np.asarray(t, np.float64)
    tc = t - t.mean()
    denom = float(tc @ tc)
    if denom <= 0.0:
        return B - B.mean(axis=1, keepdims=True)
    c = ((B - B.mean(axis=1, keepdims=True)) @ tc) / denom   # (R,)
    out = B - np.outer(c, t)
    return out - out.mean(axis=1, keepdims=True)


_FRAC_YEAR_BY_RACE: dict[int, float] | None = None


def frac_year_for_races(race_ids: np.ndarray) -> np.ndarray:
    """Fractional-year race dates for a race_id array (cached competitions lookup)."""
    global _FRAC_YEAR_BY_RACE
    if _FRAC_YEAR_BY_RACE is None:
        c = pd.read_parquet(DATA_DIR / "competitions.parquet",
                            columns=["race_id", "date"])
        _FRAC_YEAR_BY_RACE = dict(zip(c["race_id"].to_numpy(),
                                      frac_years(c["date"].to_numpy())))
    return np.array([_FRAC_YEAR_BY_RACE[int(r)] for r in race_ids], np.float64)


def poly_resid_matrix(t: np.ndarray, deg: int = 2) -> np.ndarray:
    """Projection M (n, n) that removes span{1, t_c, ..., t_c^deg}.

    Residual of a contrast vector y is `y @ M`; for a (R, n) bootstrap matrix
    of contrasts, `Y @ M` residualizes every replicate in one shot. Symmetric
    idempotent, so it applies identically on either side. Used by the cross-
    slice comparison to strip the convention-pinned smooth-in-date gauge
    family {G0, G1, Gq} from a between-slice v_j difference before ranking.
    """
    t_c = t - t.mean()
    X = np.column_stack([t_c ** k for k in range(deg + 1)])
    H = X @ np.linalg.pinv(X)
    return np.eye(len(t_c)) - H


def rank_fast(x: np.ndarray) -> np.ndarray:
    """Ascending rank, 1 = smallest value (= fastest race)."""
    return np.argsort(np.argsort(x)) + 1


def rank_matrix(M: np.ndarray) -> np.ndarray:
    """Per-row ranks of a (R, n) matrix; rank 1 = LARGEST value in the row."""
    return np.argsort(np.argsort(-M, axis=1), axis=1) + 1


def point_ranks(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(rank_desc, rank_asc): 1 = largest / 1 = smallest."""
    desc = np.argsort(np.argsort(-x)) + 1
    asc = np.argsort(np.argsort(x)) + 1
    return desc, asc


def pct(a: np.ndarray, q, axis=0) -> np.ndarray:
    return np.percentile(a, q, axis=axis)


def tier_labels(p: np.ndarray | None, point_rank: np.ndarray, n_top: int,
                p_headline: float, p_tie: float) -> np.ndarray:
    """'headline' / 'tie' / '' from P(top-N); falls back to 'point' (point
    rank <= N) when no bootstrap is available."""
    if p is None:
        return np.where(point_rank <= n_top, "point", "")
    return np.where(p >= p_headline, "headline",
                    np.where(p >= p_tie, "tie", ""))
