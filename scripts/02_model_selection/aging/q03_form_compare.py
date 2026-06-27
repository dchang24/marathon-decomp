"""Tabulate the aging-form comparison numbers -- recomputed from CSVs, no re-fit.

Pulls the key model-selection numbers into one console table per (slice, nu):

  * in-sample   -- loglik, AIC, BIC, n_params           (grid/metrics.csv, best init)
  * held-out CV -- CV log-density / cell, held-out RMSE  (cv/form_selection.csv)
  * gap vs best -- dCV/cell = (CV-best form) - (this form)
  * paired z    -- significance of that gap across folds (cv/cv_folds.csv)

and appends the **baseline (no-aging) nu=8** row from
``results/model_selection/baseline/`` so you can see how big the whole aging block
is relative to the differences *between* aging forms.

PAIRED z -- what it is
----------------------
For two forms A (the CV-best) and C that share the *same* K folds, let
``d_f = percell_A,f - percell_C,f`` be the per-fold difference in held-out
log-density per cell. Then ``z = mean_f(d_f) / (sd_f / sqrt(K))``. It is *paired*
because each fold scores an identical set of held-out cells under both models, so
fold-difficulty cancels and we test only the consistent A-vs-C gap. |z| >~ 2 means
the gap is large relative to fold-to-fold noise; tiny dCV with large z just means a
real-but-trivial difference made "significant" by the huge cell count.

NOT every form is CV-scored: e02 by default scored only poly2/poly3 + spline3..6
(gscalar/gvarying). For forms without CV (e.g. poly4/poly6, any goff) the CV
columns are blank and only the in-sample numbers are shown.

The baseline CV uses the *same* slice + fold seed + K as the aging stage (identical
per-fold n_test), so the baseline-vs-aging gap is directly comparable -- the paired
z for the baseline row is the aging block's held-out value-add.

Output: console tables + a tidy rollup CSV.
  -> results/model_selection/aging/form_compare.csv

Run::

    python scripts/02_model_selection/aging/q03_form_compare.py
    python scripts/02_model_selection/aging/q03_form_compare.py --slices ALL_M ALL_B --nu 8 inf
    python scripts/02_model_selection/aging/q03_form_compare.py --cands spline4-gvarying spline5-gvarying
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # scripts/

from marathon_decomp.config import RESULTS_DIR  # noqa: E402

OUT_ROOT = RESULTS_DIR / "model_selection" / "aging"
BASE_ROOT = RESULTS_DIR / "model_selection" / "baseline"

# forms to show by default: the flexible gvarying picks, plus spline4 in all three
# gamma flavours so the gamma axis is visible too.
DEFAULT_CANDS = [
    "poly4-gvarying", "poly6-gvarying",
    "spline4-gvarying", "spline5-gvarying", "spline6-gvarying",
    "spline4-gscalar", "spline4-goff",
]


def _f(df: pd.DataFrame, nu: float) -> pd.DataFrame:
    return df[np.isclose(df["nu"].astype(float), nu)]


def _insample(slug: str, nu: float) -> pd.DataFrame:
    """Best-init (max loglik) in-sample row per cand from grid/metrics.csv."""
    m = _f(pd.read_csv(OUT_ROOT / slug / "grid" / "metrics.csv"), nu)
    idx = m.groupby("cand")["loglik"].idxmax()
    return m.loc[idx, ["cand", "loglik", "aic", "bic", "n_params"]].set_index("cand")


def _cv(slug: str, nu: float) -> pd.DataFrame | None:
    p = OUT_ROOT / slug / "cv" / "form_selection.csv"
    if not p.is_file():
        return None
    return _f(pd.read_csv(p), nu).set_index("cand")


def _fold_pivot(slug: str, nu: float) -> pd.DataFrame | None:
    """fold x cand matrix of held-out per-cell log-density (aging stage)."""
    p = OUT_ROOT / slug / "cv" / "cv_folds.csv"
    if not p.is_file():
        return None
    f = _f(pd.read_csv(p), nu)
    return f.pivot_table(index="fold", columns="cand", values="heldout_per_cell")


def _baseline(slug: str, nu: float):
    """(in-sample+cv row dict, per-fold percell Series) for the no-aging model."""
    nsp = BASE_ROOT / slug / "nu_selection.csv"
    cfp = BASE_ROOT / slug / "cv_folds.csv"
    if not nsp.is_file():
        return None, None
    ns = _f(pd.read_csv(nsp), nu)
    if "solver" in ns:
        ns = ns[ns["solver"] == "anderson"]
    if "source" in ns and (ns["source"] == "grid").any():
        ns = ns[ns["source"] == "grid"]
    if ns.empty:
        return None, None
    r = ns.iloc[0]
    row = dict(loglik=r.get("full_loglik"), aic=r.get("full_aic"), bic=r.get("full_bic"),
               n_params=np.nan, cv_per_cell=r.get("cv_per_cell"),
               heldout_rmse=r.get("heldout_rmse"))
    fold = None
    if cfp.is_file():
        cf = _f(pd.read_csv(cfp), nu)
        if "solver" in cf and (cf["solver"] == "anderson").any():
            cf = cf[cf["solver"] == "anderson"]
        fold = cf.groupby("fold")["heldout_per_cell"].mean()
    return row, fold


def _paired_z(best_fold: pd.Series, cand_fold: pd.Series) -> float:
    """z of mean per-fold (best - cand), paired on the shared fold index."""
    if best_fold is None or cand_fold is None:
        return np.nan
    d = (best_fold - cand_fold).dropna()
    if len(d) < 2:
        return np.nan
    se = d.std(ddof=1) / np.sqrt(len(d))
    return float(d.mean() / se) if se > 0 else np.nan


def run_slice(slug: str, nu: float, cands: list[str]) -> pd.DataFrame:
    ins = _insample(slug, nu)
    cv = _cv(slug, nu)
    fp = _fold_pivot(slug, nu)
    base_row, base_fold = _baseline(slug, nu)

    # CV-best aging form (the reference for gap / z)
    best = cv["cv_per_cell"].idxmax() if cv is not None and not cv.empty else None
    best_cv = float(cv.loc[best, "cv_per_cell"]) if best else np.nan
    best_fold = fp[best] if (fp is not None and best in (fp.columns if fp is not None else [])) else None
    n_test = int(cv.loc[best, "n_test_total"]) if best else 0

    rows = []
    show = list(dict.fromkeys(cands + ([best] if best and best not in cands else [])))
    for c in show:
        rec = dict(slug=slug, nu=("inf" if not np.isfinite(nu) else f"{nu:g}"), form=c,
                   is_best=(c == best))
        if c in ins.index:
            rec.update(k=int(ins.loc[c, "n_params"]), loglik=ins.loc[c, "loglik"],
                       aic=ins.loc[c, "aic"], bic=ins.loc[c, "bic"])
        if cv is not None and c in cv.index:
            ccv = float(cv.loc[c, "cv_per_cell"])
            rec.update(cv_per_cell=ccv, heldout_rmse=float(cv.loc[c, "heldout_rmse"]),
                       dCV_vs_best=best_cv - ccv,
                       paired_z=_paired_z(best_fold, fp[c] if c in fp.columns else None))
        rows.append(rec)

    # baseline (no aging) row
    if base_row is not None:
        bcv = base_row["cv_per_cell"]
        rows.append(dict(slug=slug, nu=("inf" if not np.isfinite(nu) else f"{nu:g}"),
                         form="[baseline no-aging]", is_best=False,
                         loglik=base_row["loglik"], aic=base_row["aic"], bic=base_row["bic"],
                         cv_per_cell=bcv, heldout_rmse=base_row["heldout_rmse"],
                         dCV_vs_best=(best_cv - bcv) if np.isfinite(best_cv) else np.nan,
                         paired_z=_paired_z(best_fold, base_fold)))

    df = pd.DataFrame(rows)
    cols = ["form", "is_best", "k", "loglik", "aic", "bic",
            "cv_per_cell", "dCV_vs_best", "paired_z", "heldout_rmse"]
    df = df.reindex(columns=["slug", "nu"] + cols)

    print(f"\n=== {slug}   nu={'inf' if not np.isfinite(nu) else nu:g}   "
          f"CV-best={best}   (K-fold n_test={n_test:,}) ===")
    disp = df.drop(columns=["slug", "nu"]).copy()
    with pd.option_context("display.float_format", lambda v: f"{v:,.6g}",
                           "display.width", 200, "display.max_columns", 20):
        print(disp.to_string(index=False, na_rep="-"))
    return df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slices", nargs="+", default=["all"])
    ap.add_argument("--nu", nargs="+", default=["8"], help="nu values, e.g. 8 inf.")
    ap.add_argument("--cands", nargs="+", default=DEFAULT_CANDS)
    args = ap.parse_args()

    if args.slices == ["all"]:
        slugs = sorted(p.parent.parent.name
                       for p in OUT_ROOT.glob("*/cv/form_selection.csv"))
    else:
        all_slugs = sorted(p.parent.parent.name
                           for p in OUT_ROOT.glob("*/cv/form_selection.csv"))
        slugs = [g for s in args.slices for g in all_slugs if g == s or g.startswith(s)]
        slugs = sorted(dict.fromkeys(slugs))
    nus = [float("inf") if x.lower().startswith("inf") else float(x) for x in args.nu]

    parts = []
    for slug in slugs:
        for nu in nus:
            parts.append(run_slice(slug, nu, args.cands))
    if parts:
        out = pd.concat(parts, ignore_index=True)
        out.to_csv(OUT_ROOT / "form_compare.csv", index=False)
        print(f"\nRollup -> {OUT_ROOT / 'form_compare.csv'}")


if __name__ == "__main__":
    main()
