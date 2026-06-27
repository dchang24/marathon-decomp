"""Shared plumbing for the d_i v_j-bias validation (scripts/04_validation/d_i).

The bias test asks whether the race factor ``v_j`` absorbs the field's
career-stage composition when the model has **no** ``d_i``, and whether adding
``d_i`` removes it. The pipeline is split into independent steps so a bug in one
delta estimator never forces recomputing the other (see this dir's README):

    q01_delta_eb.py   -- producer: delta_j^EB   from the `full` fit          -> delta_eb.csv
    q02_delta_loo.py  -- producer: delta_j^LOO  from `agingS4gv` residuals    -> delta_loo.csv
    q03_vj_bias.py    -- consumer: year-partialled corr(delta, v) per fit     -> vj_bias.{md,csv}
    p01_vj_bias.py    -- plot: residualized delta-vs-v scatter

This module holds only what those scripts share: fit locators, the
within-athlete-centered career-age regressor, per-race year, the year-linear
partialling, small (scipy-free) stats, the delta-table merge I/O, and the
console+markdown ``Report`` (same house style as the aging_vs_drift study).

``delta_j`` definitions (both averaged over the finishers of race j; an
ineligible athlete contributes 0, matching how the model treats it):

  * EB  : delta_j = mean_{i in j} ( d_i^full * an_tilde_ij )         [uses `full`]
  * LOO : delta_j = mean_{i in j} ( b_i^{(-j)} * an_tilde_ij )       [uses `no-d` resid]
          where b_i is athlete i's OLS slope of the no-d residual on
          an_tilde, and b_i^{(-j)} leaves race j's own residual out (rank-1
          downdate) so delta_j shares no information with v_j^{no-d}.

``an_tilde_ij = A_n_ij - mean_i(A_n)`` is the model's drift regressor
(``Model.An_tilde_cell``); it is recomputed here directly from ``fd`` so the
producers do not depend on a model's internals for the regressor itself.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # scripts/

from marathon_decomp import ModelConfig, registry                # noqa: E402
from marathon_decomp.config import RESULTS_DIR                    # noqa: E402

from aging_common.fitting import aging_cfg, aging_stem            # noqa: E402
from drift_common.fitting import drift_cfg, drift_stem            # noqa: E402

OUT_ROOT = RESULTS_DIR / "validation" / "d_i"
MODELS_ROOT = RESULTS_DIR / "models"

# delta-table schema (one row per race); slug already encodes (slice, mrc).
DELTA_COLS = ["slug", "slice", "mrc", "nu", "race_idx", "race_id", "year",
              "delta", "n_finishers", "n_eligible_in_race"]


# -- locating the two registered fits the test consumes -------------------
def aging_dir(spec, nu: float) -> Path:
    """The slice's registered no-d fit (agingS4gv_nu8p00_best)."""
    cfg = aging_cfg(nu)
    parent = MODELS_ROOT / registry.slice_slug(spec)
    return registry.fit_path(parent, aging_stem(cfg), spec, cfg, resample_tag="base")


def full_dir(spec, nu: float) -> Path:
    """The slice's registered +d fit (full_nu8p00_best)."""
    cfg = drift_cfg(nu, variant="full")
    parent = MODELS_ROOT / registry.slice_slug(spec)
    return registry.fit_path(parent, drift_stem("full", nu), spec, cfg, resample_tag="base")


def present(fit_dir: Path) -> bool:
    return (fit_dir / "fit.pkl").is_file() and (fit_dir / "manifest.json").is_file()


# -- per-obs / per-race quantities straight off the FitData ---------------
def an_tilde(fd) -> np.ndarray:
    """Within-athlete-centered career age (length N) = A_n - mean_i(A_n)."""
    row = fd.row_idx
    n = np.bincount(row, minlength=fd.I)
    s = np.bincount(row, weights=fd.A_n, minlength=fd.I)
    bar = np.where(n > 0, s / np.maximum(n, 1.0), 0.0)
    return fd.A_n - bar[row]


def n_per_athlete(fd) -> np.ndarray:
    return np.bincount(fd.row_idx, minlength=fd.I)


def race_year(fd) -> np.ndarray:
    """Per-race calendar year (length J), float for the partialling design."""
    return pd.DatetimeIndex(fd.race_date).year.to_numpy().astype(float)


def per_race_mean(values: np.ndarray, col: np.ndarray, J: int) -> np.ndarray:
    """Mean of a per-obs quantity within each race (0 where the race is empty)."""
    tot = np.bincount(col, weights=values, minlength=J)
    n = np.bincount(col, minlength=J)
    return np.where(n > 0, tot / np.maximum(n, 1.0), 0.0)


def yhat(model) -> np.ndarray:
    """Full per-obs prediction = sum of the model's decomposition terms."""
    return sum(model._yhat_terms().values())


_TINY = 1e-12


def loo_delta(fd, r: np.ndarray, an: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Leave-one-out LOO delta_j from a per-obs residual `r` and regressor `an`.

    Per athlete, OLS slope of r on an (an is within-athlete-centered, mean 0);
    the leave-one-out slope drops the obs's own residual via rank-1 downdate so
    delta_j shares no information with that race's residual. Returns
    (delta[J], eligible_athlete[I]). Eligible = n_i>=2 and sum(an^2)>0.

    This is the single source of truth for q02 (observed) and q04 (permuted r).
    """
    row, col = fd.row_idx, fd.col_idx
    S_ar = np.bincount(row, weights=an * r, minlength=fd.I)
    S_aa = np.bincount(row, weights=an * an, minlength=fd.I)
    n_i = np.bincount(row, minlength=fd.I)
    elig = (n_i >= 2) & (S_aa > _TINY)

    num = S_ar[row] - an * r
    den = S_aa[row] - an * an
    ok = elig[row] & (den > _TINY)
    b_loo = np.where(ok, num / np.where(den > _TINY, den, 1.0), 0.0)
    contrib = b_loo * an
    return per_race_mean(contrib, col, fd.J), elig


def permute_within_athlete(r: np.ndarray, row: np.ndarray,
                           rng: np.random.Generator) -> np.ndarray:
    """Shuffle r within each athlete (row group), preserving group membership.

    Breaks the within-athlete an<->r pairing (the trajectory) while keeping each
    athlete's residual multiset and the race assignment intact -- the null for
    the v_j-bias permutation test.
    """
    keys = rng.random(r.size)
    shuffled_slots = np.lexsort((keys, row))   # positions grouped by row, random within
    grouped = np.argsort(row, kind="stable")   # positions grouped by row, original order
    out = np.empty_like(r)
    out[shuffled_slots] = r[grouped]
    return out


# -- year-linear partialling + scipy-free stats ---------------------------
def partial_on_year(x: np.ndarray, year: np.ndarray) -> np.ndarray:
    """Residual of x after OLS regression on [1, year] (date control)."""
    x = np.asarray(x, float)
    Z = np.column_stack([np.ones(year.size), np.asarray(year, float)])
    beta, *_ = np.linalg.lstsq(Z, x, rcond=None)
    return x - Z @ beta


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, float); b = np.asarray(b, float)
    if a.size < 3 or a.std() < 1e-12 or b.std() < 1e-12:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def _rank(a: np.ndarray) -> np.ndarray:
    return np.asarray(a, float).argsort().argsort().astype(float)


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    return pearson(_rank(a), _rank(b))


def ols_slope(y: np.ndarray, x: np.ndarray) -> float:
    """Slope of y ~ x (both already mean-zero residuals); cov(x,y)/var(x)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    vx = float(np.var(x))
    return float(np.cov(x, y)[0, 1] / vx) if vx > 1e-18 else float("nan")


# -- delta-table merge I/O (rerun one slice without losing the others) ----
def delta_path(estimator: str) -> Path:
    return OUT_ROOT / f"delta_{estimator}.csv"


def merge_delta(new_df: pd.DataFrame, estimator: str) -> Path:
    """Upsert new_df into delta_<estimator>.csv, keyed on slug (mrc-aware).

    Rows for any slug present in new_df are replaced; other slugs are kept, so
    running one slice then another accumulates rather than overwrites.
    """
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    new_df = new_df[DELTA_COLS]
    p = delta_path(estimator)
    if p.is_file():
        old = pd.read_csv(p)
        old = old[~old["slug"].isin(new_df["slug"].unique())]
        out = pd.concat([old, new_df], ignore_index=True)
    else:
        out = new_df
    out = out.sort_values(["slug", "race_idx"]).reset_index(drop=True)
    out.to_csv(p, index=False)
    return p


def read_delta(estimator: str) -> pd.DataFrame | None:
    p = delta_path(estimator)
    return pd.read_csv(p) if p.is_file() else None


def upsert_csv(path: Path, new_df: pd.DataFrame, keys: list[str]) -> Path:
    """Insert/replace new_df rows into path, keyed on `keys` (durable accumulator).

    Rows in the existing file whose key tuple appears in new_df are replaced;
    all others are kept. Used by q03/q04 so a per-cell or single-slice run never
    wipes the other slices already in the one-file-per-artifact CSV.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file():
        old = pd.read_csv(path)
        kdf = new_df[keys].drop_duplicates()
        marked = old.merge(kdf, on=keys, how="left", indicator=True)
        old = old[marked["_merge"].to_numpy() == "left_only"]
        out = pd.concat([old, new_df], ignore_index=True)
    else:
        out = new_df
    out.to_csv(path, index=False)
    return path


# -- console + markdown report (same style as aging_vs_drift) -------------
def _md_cell(x) -> str:
    return str(x).replace("|", "\\|")


def _df_to_md(df: pd.DataFrame, *, index: bool) -> str:
    cols = ([df.index.name or ""] if index else []) + [str(c) for c in df.columns]
    rows = []
    for idx, row in df.iterrows():
        cells = ([_md_cell(idx)] if index else []) + [_md_cell(v) for v in row]
        rows.append("| " + " | ".join(cells) + " |")
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    return "\n".join([header, sep] + rows)


class Report:
    """Echoes everything to the console and accumulates a markdown mirror."""

    def __init__(self) -> None:
        self.md: list[str] = []

    def line(self, text: str = "") -> None:
        print(text, flush=True)
        self.md.append(text)

    def head(self, text: str, level: int = 2) -> None:
        print(f"\n{'=' * 4} {text} {'=' * 4}" if level == 2 else f"\n-- {text} --",
              flush=True)
        self.md.append("")
        self.md.append(f"{'#' * level} {text}")
        self.md.append("")

    def table(self, df: pd.DataFrame, *, index: bool = False,
              floatfmt: str = "{:.6g}") -> None:
        def fmt(x):
            if isinstance(x, float):
                return "-" if not np.isfinite(x) else floatfmt.format(x)
            return str(x)
        disp = df.copy()
        for c in disp.columns:
            disp[c] = disp[c].map(fmt)
        with pd.option_context("display.width", 200, "display.max_columns", 40):
            print(disp.to_string(index=index), flush=True)
        self.md.append(_df_to_md(disp, index=index))
        self.md.append("")

    def save(self, path: Path, title: str) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("# " + title + "\n\n" + "\n".join(self.md) + "\n",
                        encoding="utf-8")
        print(f"\nmarkdown -> {path}")
