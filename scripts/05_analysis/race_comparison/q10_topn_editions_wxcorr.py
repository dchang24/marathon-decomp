"""Top-N slowest / fastest race editions by the WEATHER-CORRECTED v_j.

Net-of-weather counterpart to q02 (the variant flagged "deferred to the
covariate analysis" in the README). Instead of v_j -- which bundles course +
that day's weather + field -- it ranks ``v_corr_j`` from
covariate/q07_weather_corrected_vj.py::

    v_corr_j = v_j - b_w * (temp_field_j - temp_ref)

i.e. each edition normalized to a reference air temperature using the partial
air-temp slope ``b_w`` of the joint v_j ~ temp + total_gain regression. ONLY
air temperature is removed; course and field stay in. So the estimand is:

    "slowest / fastest RELATIVE TO ITS ERA, *at reference air temperature*,
     still bundling course + field + the unexplained part."

Reading vs q02: a race that is hot-edition-driven in q02 should move toward the
middle here; a course that is genuinely hard at neutral weather should hold its
rank. The contrast between the two tables is the whole point.

SELECTION BY RANK STABILITY (same as q02). For every q07 bootstrap replicate
(which already propagates BOTH v_j sampling noise and b_w uncertainty -- q07
refits b_w per replicate) we rank and accumulate P(top-N):

    headline  P >= p_headline (default 0.50)
    tie       p_tie <= P < p_headline (default 0.25)

With no bootstrap for a slice, falls back to the point top-N (tier 'point').
The --min-n floor drops races with fewer than that many finishers in the slice
before ranking (default 0). Only weather-COVERED editions (air temp present)
are correctable and hence enter this ranking; coverage vs q02 is reported.

Outputs, one dir per slice slug (results/analysis/race_comparison/{slug}/),
filenames stamped with parameters so runs never overwrite each other::

    topn_editions_wxcorr_minn{min_n}_top{n_top}.csv   all covered races
    topn_editions_wxcorr_minn{min_n}_top{n_top}.md    headline + tie tables
    + the same tables to stdout.

Run after covariate/q07_weather_corrected_vj.py ::

    python scripts/05_analysis/race_comparison/q04_topn_editions_wxcorr.py
    python scripts/05_analysis/race_comparison/q04_topn_editions_wxcorr.py \
        --slices ALL_B Po10_M --min-n 25 --n-top 10
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from marathon_decomp.config import DATA_DIR, RESULTS_DIR  # noqa: E402
from race_common import (  # noqa: E402
    OUT_ROOT, v_to_min, pct, point_ranks, rank_matrix, tier_labels,
)
from q02_topn_editions import report_slice  # noqa: E402  (reuse table formatting)

COVARIATE_DIR = RESULTS_DIR / "analysis" / "covariate" / "07_weather_corrected_vj"
DEFAULT_SLICES = ["ALL_B", "ALL_M", "ALL_W", "Po10_B", "Po10_M", "Po10_W"]


def slice_slug(slice_name: str, window: str, mrc: str) -> str:
    """Covariate label 'ALL_B' -> q02-style slug 'ALL_B_14-25_mrc2'."""
    return f"{slice_name}_{window}_{mrc}"


def boot_matrix(sub: pd.DataFrame, race_ids: np.ndarray,
                point_v: np.ndarray) -> np.ndarray | None:
    """q07 corrected replicates (already one slice) -> (R, n) aligned, or None."""
    if sub.empty:
        return None
    wide = (sub.pivot(index="run_id", columns="race_id", values="v_corr")
            .reindex(columns=race_ids))
    if wide.empty:
        return None
    B = wide.to_numpy(np.float64)
    bad = ~np.isfinite(B)
    if bad.any():
        print(f"    [warn] {int(bad.any(axis=0).sum())} races missing from "
              f"corrected bootstrap; filled with point v_corr")
        B = np.where(bad, np.broadcast_to(point_v, B.shape), B)
    return B


def analyse_slice(slice_name: str, point: pd.DataFrame,
                  boot_all: pd.DataFrame | None, comps: pd.DataFrame, *,
                  slug: str, min_n: int, n_top: int,
                  p_headline: float, p_tie: float) -> pd.DataFrame | None:
    d = point[point["slice"] == slice_name].copy()
    if d.empty:
        print(f"  [skip] {slice_name}: no corrected v_j rows")
        return None
    d = d.drop(columns="slice")             # q07 label; replaced by the slug below
    d = d.merge(comps, on="race_id", how="left")
    d["series_key"] = d["series_key"].fillna("unknown").astype(str)
    d["year"] = d["year"].fillna(0).astype(int)

    n_cov = len(d)
    d = d[d["n_j"] >= min_n].reset_index(drop=True)
    if n_cov - len(d):
        print(f"  {slice_name}: excluded {n_cov - len(d)}/{n_cov} races with n_j < {min_n}")
    if len(d) <= n_top:
        print(f"  [skip] {slice_name}: only {len(d)} eligible races (<= n_top)")
        return None

    v = d["v_corr"].to_numpy(np.float64)
    rank_slow, rank_fast = point_ranks(v)
    d["v"] = v                                  # report_slice prints column 'v'
    d["rank_slow"], d["rank_fast"] = rank_slow, rank_fast
    d["mult"] = np.exp(v)
    d["min_at_3h"] = v_to_min(v)

    B = None
    if boot_all is not None:
        bs = boot_all[boot_all["slice"] == slice_name]
        if not bs.empty:
            B = boot_matrix(bs, d["race_id"].to_numpy(np.int64), v)

    p_slow = p_fast = None
    if B is not None:
        d["v_lo95"], d["v_hi95"] = pct(B, 2.5), pct(B, 97.5)
        rs, rf = rank_matrix(B), rank_matrix(-B)
        p_slow, p_fast = (rs <= n_top).mean(axis=0), (rf <= n_top).mean(axis=0)
        d["p_slow"], d["p_fast"] = p_slow, p_fast
        d["rank_slow_lo"], d["rank_slow_hi"] = pct(rs, 2.5), pct(rs, 97.5)
        d["rank_fast_lo"], d["rank_fast_hi"] = pct(rf, 2.5), pct(rf, 97.5)
        d["n_boot"] = B.shape[0]
    else:
        print(f"    [note] no corrected bootstrap for {slice_name}; point-only tiers")
        for c in ("v_lo95", "v_hi95", "p_slow", "p_fast", "rank_slow_lo",
                  "rank_slow_hi", "rank_fast_lo", "rank_fast_hi"):
            d[c] = np.nan
        d["n_boot"] = 0

    d["tier_slow"] = tier_labels(p_slow, rank_slow, n_top, p_headline, p_tie)
    d["tier_fast"] = tier_labels(p_fast, rank_fast, n_top, p_headline, p_tie)

    d.insert(0, "slice", slug)
    d.insert(1, "model", "full_wxcorr")
    d.insert(2, "min_n", min_n)
    d.insert(3, "n_top", n_top)
    return d


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--slices", nargs="+", default=DEFAULT_SLICES,
                    help="covariate slice labels (ALL_B, Po10_M, ...).")
    ap.add_argument("--window", default="14-25")
    ap.add_argument("--mrc", default="mrc2")
    ap.add_argument("--n-top", type=int, default=10)
    ap.add_argument("--min-n", type=int, default=0)
    ap.add_argument("--p-headline", type=float, default=0.50)
    ap.add_argument("--p-tie", type=float, default=0.25)
    args = ap.parse_args()

    pp = COVARIATE_DIR / "weather_corrected_vj__6slices.parquet"
    if not pp.is_file():
        raise SystemExit(f"missing {pp}; run covariate/q07_weather_corrected_vj.py first")
    point = pd.read_parquet(pp)
    bp = COVARIATE_DIR / "weather_corrected_vj_boot__6slices.parquet"
    boot_all = pd.read_parquet(bp) if bp.is_file() else None
    comps = pd.read_parquet(DATA_DIR / "competitions.parquet",
                            columns=["race_id", "series_key", "country", "date", "year"])

    stem = f"topn_editions_wxcorr_minn{args.min_n}_top{args.n_top}"
    out_all, n_done = [], 0
    for s in args.slices:
        slug = slice_slug(s, args.window, args.mrc)
        print(f"analysing {slug} (weather-corrected) ...")
        d = analyse_slice(s, point, boot_all, comps, slug=slug,
                          min_n=args.min_n, n_top=args.n_top,
                          p_headline=args.p_headline, p_tie=args.p_tie)
        if d is None:
            continue
        md, out = report_slice(d, args.n_top)
        out_all += out

        out_dir = OUT_ROOT / slug
        out_dir.mkdir(parents=True, exist_ok=True)
        num = d.select_dtypes(include="number").columns
        d[num] = d[num].round(6)
        d.to_csv(out_dir / f"{stem}.csv", index=False)
        head = [
            "# Top-N slowest / fastest race editions by WEATHER-CORRECTED v_j",
            "",
            "Estimand: v_corr_j = v_j - b_w*(temp - temp_ref) -- each edition",
            "normalized to a reference air temperature (air temp only; course +",
            "field stay in). 'Slowest/fastest' = relative to era AT REFERENCE",
            "WEATHER. Selection by bootstrap rank stability P(top-N) (headline",
            "P>=0.5, tie P>=0.25); the q07 bootstrap propagates both v_j and b_w",
            "noise. min@3:00 = 180*(exp(v_corr)-1) = minutes vs the average race at 3:00:00.",
            "", f"slice: {slug}; min_n = {args.min_n}; n_top = {args.n_top}; "
            "model = full_wxcorr (air-temp corrected)", ""]
        (out_dir / f"{stem}.md").write_text("\n".join(head + md), encoding="utf-8")
        print(f"  wrote {out_dir / (stem + '.csv')}  ({len(d)} rows)")
        print(f"  wrote {out_dir / (stem + '.md')}")
        n_done += 1

    if not n_done:
        raise SystemExit("no slices produced a result")

    print("\n".join(out_all))
    print()


if __name__ == "__main__":
    main()
