"""Consumer: does d_i remove the field-composition bias in v_j?

For each delta estimator on disk (``delta_eb.csv`` from q01, ``delta_loo.csv``
from q02) and each (slice, mrc), this correlates the per-race field career-stage
tilt ``delta_j`` against the race factor ``v_j`` from two fits that differ only in
the d_i block (aging is on in both):

    no-d : v_j from agingS4gv_nu8p00_best
    +d   : v_j from full_nu8p00_best

Each correlation is reported raw and **year-partialled** (v and delta each
residualized on [1, year]; the date control). The bias test predicts:

    corr(delta, v_no-d) > 0     (the composition leaks into v when d is absent)
    corr(delta, v_+d)   ~ 0     (d_i absorbs it)

so ``drop = corr_no-d - corr_+d`` is large and ``frac_removed -> 1``. This is a
cheap read of the delta tables + the two v-vectors -- it never recomputes a
delta, so re-running after a producer fix is safe and fast.

This script does NOT establish significance; the within-athlete permutation null
(a later q04) is the proper test of "is corr_no-d real". The PASS flag here is a
descriptive heuristic only.

Outputs (under results/validation/d_i/):
  * vj_bias.md          -- the full console run, formatted
  * bias_corr.csv       -- per (slug, mrc, estimator, fit): raw/partial pearson+spearman, slope, n_races
  * bias_summary.csv    -- per (slug, mrc, estimator): corr_nod, corr_d, drop, frac_removed, pass
  * delta_agreement.csv -- per (slug, mrc): corr(delta_eb, delta_loo) where both exist

Run::

    python scripts/04_validation/d_i/q03_vj_bias.py --slices Po10_W
    python scripts/04_validation/d_i/q03_vj_bias.py                 # all 8, mrc 2 & 5
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # scripts/

from marathon_decomp import load_slice, registry                # noqa: E402
from marathon_decomp.config import display_path                  # noqa: E402

from baseline_common import slices as S                          # noqa: E402
import delta_common as DC                                        # noqa: E402  (sibling)

ESTIMATORS = ("eb", "loo")
# descriptive PASS heuristic (NOT inference; see q04 permutation null)
_POS = 0.05          # no-d partial corr must clear this to count as "positive"
_SHRINK = 0.5        # +d partial corr must be below this fraction of no-d's


def _delta_for_slug(table: pd.DataFrame | None, slug: str) -> pd.DataFrame | None:
    if table is None:
        return None
    sub = table[table["slug"] == slug]
    return sub.sort_values("race_idx") if len(sub) else None


def _corr_block(delta: np.ndarray, v: np.ndarray, year: np.ndarray) -> dict:
    dv = DC.partial_on_year(delta, year)
    vv = DC.partial_on_year(v, year)
    return dict(
        pearson_raw=DC.pearson(delta, v),
        spearman_raw=DC.spearman(delta, v),
        pearson_partial=DC.pearson(dv, vv),
        spearman_partial=DC.spearman(dv, vv),
        slope_partial=DC.ols_slope(vv, dv),
    )


def run_cell(name: str, mrc: int, nu: float, tables: dict, rep: DC.Report,
             corr_rows: list, summ_rows: list, agree_rows: list) -> bool:
    spec = S.build_spec(name, min_race_count=mrc)
    slug = registry.slice_slug(spec)
    adir, fdir = DC.aging_dir(spec, nu), DC.full_dir(spec, nu)
    if not (DC.present(adir) and DC.present(fdir)):
        return False
    have = {e: _delta_for_slug(tables.get(e), slug) for e in ESTIMATORS}
    if not any(d is not None for d in have.values()):
        return False

    fd = load_slice(spec)
    v_nod = np.asarray(registry.load_fit(adir, fd).params["v"], float)
    v_d = np.asarray(registry.load_fit(fdir, fd).params["v"], float)
    year = DC.race_year(fd)

    rep.head(f"{slug}   nu={nu:g}", level=2)
    rep.line(f"J={fd.J:,} races   v: no-d=agingS4gv, +d=full   "
             f"estimators present: {[e for e in ESTIMATORS if have[e] is not None]}")

    # delta_eb vs delta_loo agreement (when both present)
    if have["eb"] is not None and have["loo"] is not None:
        de = have["eb"].set_index("race_idx")["delta"]
        dl = have["loo"].set_index("race_idx")["delta"]
        j = de.index.intersection(dl.index)
        ag = DC.pearson(de.loc[j].to_numpy(), dl.loc[j].to_numpy())
        agree_rows.append(dict(slug=slug, slice=name, mrc=mrc, nu=float(nu),
                               corr_eb_loo=ag, sd_eb=float(de.std()), sd_loo=float(dl.std())))
        rep.line(f"\ndelta agreement  corr(delta_EB, delta_LOO) = {ag:.4f}")

    for est in ESTIMATORS:
        dsub = have[est]
        if dsub is None:
            continue
        idx = dsub["race_idx"].to_numpy()
        delta = dsub["delta"].to_numpy()
        bn = _corr_block(delta, v_nod[idx], year[idx])
        bd = _corr_block(delta, v_d[idx], year[idx])

        rep.head(f"estimator = {est.upper()}", level=3)
        tbl = pd.DataFrame({
            "quantity": ["pearson_raw", "spearman_raw", "pearson_partial",
                         "spearman_partial", "slope_partial"],
            "no-d (agingS4gv)": [bn[k] for k in ("pearson_raw", "spearman_raw",
                                                 "pearson_partial", "spearman_partial", "slope_partial")],
            "+d (full)": [bd[k] for k in ("pearson_raw", "spearman_raw",
                                          "pearson_partial", "spearman_partial", "slope_partial")],
        })
        rep.table(tbl, floatfmt="{:.4f}")

        cn, cd = bn["pearson_partial"], bd["pearson_partial"]
        drop = cn - cd
        frac = float(1.0 - cd / cn) if np.isfinite(cn) and abs(cn) > 1e-9 else float("nan")
        passed = bool(np.isfinite(cn) and cn > _POS and abs(cd) < _SHRINK * abs(cn))
        rep.line(f"\nyear-partialled pearson:  no-d={cn:+.4f}  +d={cd:+.4f}  "
                 f"drop={drop:+.4f}  frac_removed={frac:.3f}  -> "
                 f"{'PASS' if passed else 'check'}")

        for fit, blk in (("no-d", bn), ("+d", bd)):
            corr_rows.append(dict(slug=slug, slice=name, mrc=mrc, nu=float(nu),
                                  estimator=est, fit=fit, n_races=int(fd.J), **blk))
        summ_rows.append(dict(slug=slug, slice=name, mrc=mrc, nu=float(nu),
                              estimator=est, corr_nod=cn, corr_d=cd, drop=drop,
                              frac_removed=frac, pass_heuristic=passed))
    return True


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--slices", nargs="+", default=["all"])
    ap.add_argument("--mrc", type=int, nargs="+", default=[2, 5])
    ap.add_argument("--nu", type=float, default=8.0)
    args = ap.parse_args()

    names = S.resolve_names(args.slices, ap)
    mrcs = sorted(dict.fromkeys(args.mrc))
    tables = {e: DC.read_delta(e) for e in ESTIMATORS}
    if all(t is None for t in tables.values()):
        ap.error("no delta tables found; run q01_delta_eb.py and/or q02_delta_loo.py first.")

    rep = DC.Report()
    rep.line(f"d_i v_j-bias test   nu={args.nu:g}   mrc={mrcs}   "
             f"slices: {', '.join(names)}")
    rep.line(f"delta tables present: {[e for e in ESTIMATORS if tables[e] is not None]}")
    rep.line("test: corr(delta, v_no-d) > 0  and  corr(delta, v_+d) ~ 0 "
             "(year-partialled). Significance deferred to q04 permutation null.")

    corr_rows, summ_rows, agree_rows = [], [], []
    n_cells = 0
    for name in names:
        for mrc in mrcs:
            if run_cell(name, mrc, args.nu, tables, rep, corr_rows, summ_rows, agree_rows):
                n_cells += 1

    if n_cells == 0:
        rep.line("\nNo (slice, mrc) had both fits + a delta table. Run the "
                 "producers and the e01_fit_drift.py fits first.")
        rep.save(DC.OUT_ROOT / "vj_bias.md", "d_i v_j-bias test")
        return

    # cross-slice headline
    if summ_rows:
        rep.head("Headline: year-partialled pearson across slices", level=2)
        sdf = pd.DataFrame(summ_rows)
        rep.table(sdf[["slug", "estimator", "corr_nod", "corr_d", "drop",
                       "frac_removed", "pass_heuristic"]], floatfmt="{:.4f}")

    DC.OUT_ROOT.mkdir(parents=True, exist_ok=True)
    # upsert keyed on slug so a single-slice run never wipes the other slices
    DC.upsert_csv(DC.OUT_ROOT / "bias_corr.csv", pd.DataFrame(corr_rows), keys=["slug"])
    DC.upsert_csv(DC.OUT_ROOT / "bias_summary.csv", pd.DataFrame(summ_rows), keys=["slug"])
    if agree_rows:
        DC.upsert_csv(DC.OUT_ROOT / "delta_agreement.csv", pd.DataFrame(agree_rows), keys=["slug"])
    rep.line(f"\nCSV -> {display_path(DC.OUT_ROOT)}")
    rep.line("  bias_corr.csv  bias_summary.csv"
             + ("  delta_agreement.csv" if agree_rows else ""))
    rep.save(DC.OUT_ROOT / "vj_bias.md", "d_i v_j-bias test")
    print(f"\nDone: {n_cells} (slice, mrc) cell(s) tested.")


if __name__ == "__main__":
    main()
