"""Top-N fastest / slowest race SERIES by median WEATHER-CORRECTED edition v_j.

Net-of-weather counterpart to q03. q03's "PB-chaser index" deliberately INCLUDES
each series' typical weather (London's reliably cool April is legitimately part
of why it runs fast). This script instead evaluates every series at a COMMON
reference air temperature, using the weather-corrected edition factor
``v_corr_j = v_j - b_w*(temp_field_j - temp_ref)`` from
covariate/q07_weather_corrected_vj.py (air temp only; course + field stay in):

    m_s = median_j( v_corr_j ),     over covered editions with n_j >= min_n.

ESTIMAND. "If you caught a reference-temperature (good-weather) day, how fast
does this series run relative to its era, net of course + field?" A series that
is fast only because it is typically cold loses that credit (its median is
pulled UP/slower); a series that is typically hot is no longer penalized for it
(median pulled DOWN/faster). This is the index a PB-chaser wants when picking a
race and hoping for good conditions -- the inverse of q03's "typical weather
included" reading.

LEVEL vs ORDER. temp_ref (q07's per-slice mean air temp) only shifts every
series' median by the same constant, so the RANKING is temp_ref-invariant; only
the min@3:00 *level* depends on it. Because the correction is mean-zero over the
covered editions, the corrected scale stays centred like q03's ("minutes vs the
average race"), so the two tables' min@3:00 columns are directly comparable.

Everything else matches q03: the median discards one freak edition with no
hand-picked exclusions (hence --min-editions, default 3); bootstrap medians are
recomputed INSIDE each q07 replicate (which already carries both v_j and b_w
noise), preserving the within-replicate edition correlation; selection is by
rank stability P(top-N) (headline >= 0.50, tie >= 0.25); no total order claimed.
Only weather-COVERED editions (air temp present) enter -- coverage is reported.

Outputs, one dir per slice slug (results/analysis/race_comparison/{slug}/):
    series_ranking_wxcorr_minn{min_n}_mined{min_editions}_top{n_top}.{csv,md}
    + the same tables to stdout.

Run after covariate/q07_weather_corrected_vj.py ::

    python scripts/05_analysis/race_comparison/q05_series_ranking_wxcorr.py
    python scripts/05_analysis/race_comparison/q05_series_ranking_wxcorr.py \
        --slices ALL_B Po10_M --min-n 25 --min-editions 3 --n-top 20
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
from q03_series_ranking import report_slice  # noqa: E402  (reuse table formatting)

COVARIATE_DIR = RESULTS_DIR / "analysis" / "covariate" / "07_weather_corrected_vj"
DEFAULT_SLICES = ["ALL_B", "ALL_M", "ALL_W", "Po10_B", "Po10_M", "Po10_W"]


def slice_slug(slice_name: str, window: str, mrc: str) -> str:
    return f"{slice_name}_{window}_{mrc}"


def boot_matrix(sub: pd.DataFrame, race_ids: np.ndarray,
                point_v: np.ndarray) -> np.ndarray | None:
    """q07 corrected replicates (one slice) -> (R, n) aligned to race_ids."""
    if sub.empty:
        return None
    wide = (sub.pivot(index="run_id", columns="race_id", values="v_corr")
            .reindex(columns=race_ids))
    if wide.empty:
        return None
    B = wide.to_numpy(np.float64)
    bad = ~np.isfinite(B)
    if bad.any():
        print(f"    [warn] {int(bad.any(axis=0).sum())} editions missing from "
              f"corrected bootstrap; filled with point v_corr")
        B = np.where(bad, np.broadcast_to(point_v, B.shape), B)
    return B


def analyse_slice(slice_name: str, point: pd.DataFrame,
                  boot_all: pd.DataFrame | None, comps: pd.DataFrame, *,
                  slug: str, min_n: int, min_editions: int, n_top: int,
                  p_headline: float, p_tie: float) -> pd.DataFrame | None:
    ed = point[point["slice"] == slice_name].copy()
    if ed.empty:
        print(f"  [skip] {slice_name}: no corrected v_j rows")
        return None
    ed = ed.merge(comps, on="race_id", how="left")
    ed["series_key"] = ed["series_key"].fillna("unknown").astype(str)
    ed["year"] = ed["year"].fillna(0).astype(int)

    n_all = len(ed)
    ed = ed[ed["n_j"] >= min_n].reset_index(drop=True)
    if n_all - len(ed):
        print(f"  {slug}: excluded {n_all - len(ed)}/{n_all} editions with n_j < {min_n}")

    sizes = ed.groupby("series_key").size()
    keep = sizes[sizes >= min_editions].index
    ed_k = ed[ed["series_key"].isin(keep)].reset_index(drop=True)
    print(f"  {slug}: {len(keep)} series with >= {min_editions} covered editions "
          f"({ed['series_key'].nunique() - len(keep)} dropped)")
    if len(keep) <= n_top:
        print(f"  [skip] {slug}: only {len(keep)} eligible series (<= n_top)")
        return None

    # point summary per series (median of weather-corrected v_corr)
    g = ed_k.groupby("series_key")
    sr = pd.DataFrame({
        "series_key": g.size().index,
        "country": g["country"].first().to_numpy(),
        "k_editions": g.size().to_numpy(int),
        "yr_min": g["year"].min().to_numpy(int),
        "yr_max": g["year"].max().to_numpy(int),
        "n_j_total": g["n_j"].sum().to_numpy(int),
        "median_v": g["v_corr"].median().to_numpy(np.float64),
        "iqr_v": (g["v_corr"].quantile(0.75) - g["v_corr"].quantile(0.25)).to_numpy(np.float64),
    }).reset_index(drop=True)

    m = sr["median_v"].to_numpy(np.float64)
    sr["min_at_3h"] = v_to_min(m)
    rank_slow, rank_fast = point_ranks(m)
    sr["rank_slow"], sr["rank_fast"] = rank_slow, rank_fast

    B = None
    if boot_all is not None:
        bs = boot_all[boot_all["slice"] == slice_name]
        if not bs.empty:
            B = boot_matrix(bs, ed_k["race_id"].to_numpy(np.int64),
                            ed_k["v_corr"].to_numpy(np.float64))
    p_slow = p_fast = None
    if B is not None:
        cols = {s: np.flatnonzero((ed_k["series_key"] == s).to_numpy())
                for s in sr["series_key"]}
        M = np.column_stack([np.median(B[:, cols[s]], axis=1)
                             for s in sr["series_key"]])     # (R, S)
        sr["m_lo95"], sr["m_hi95"] = pct(M, 2.5), pct(M, 97.5)
        rs, rf = rank_matrix(M), rank_matrix(-M)
        p_slow, p_fast = (rs <= n_top).mean(axis=0), (rf <= n_top).mean(axis=0)
        sr["p_slow"], sr["p_fast"] = p_slow, p_fast
        sr["rank_slow_lo"], sr["rank_slow_hi"] = pct(rs, 2.5), pct(rs, 97.5)
        sr["rank_fast_lo"], sr["rank_fast_hi"] = pct(rf, 2.5), pct(rf, 97.5)
        sr["n_boot"] = M.shape[0]
    else:
        print(f"    [note] no corrected bootstrap for {slice_name}; point-only tiers")
        for c in ("m_lo95", "m_hi95", "p_slow", "p_fast", "rank_slow_lo",
                  "rank_slow_hi", "rank_fast_lo", "rank_fast_hi"):
            sr[c] = np.nan
        sr["n_boot"] = 0

    sr["tier_slow"] = tier_labels(p_slow, rank_slow, n_top, p_headline, p_tie)
    sr["tier_fast"] = tier_labels(p_fast, rank_fast, n_top, p_headline, p_tie)

    sr.insert(0, "slice", slug)
    sr.insert(1, "model", "full_wxcorr")
    sr.insert(2, "min_n", min_n)
    sr.insert(3, "min_editions", min_editions)
    sr.insert(4, "n_top", n_top)
    return sr


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--slices", nargs="+", default=DEFAULT_SLICES,
                    help="covariate slice labels (ALL_B, Po10_M, ...).")
    ap.add_argument("--window", default="14-25")
    ap.add_argument("--mrc", default="mrc2")
    ap.add_argument("--n-top", type=int, default=20)
    ap.add_argument("--min-n", type=int, default=0)
    ap.add_argument("--min-editions", type=int, default=3)
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

    stem = (f"series_ranking_wxcorr_minn{args.min_n}"
            f"_mined{args.min_editions}_top{args.n_top}")
    out_all, n_done = [], 0
    for s in args.slices:
        slug = slice_slug(s, args.window, args.mrc)
        print(f"analysing {slug} (weather-corrected series) ...")
        sr = analyse_slice(s, point, boot_all, comps, slug=slug,
                           min_n=args.min_n, min_editions=args.min_editions,
                           n_top=args.n_top, p_headline=args.p_headline,
                           p_tie=args.p_tie)
        if sr is None:
            continue
        md, out = report_slice(sr, args.n_top)
        out_all += out

        out_dir = OUT_ROOT / slug
        out_dir.mkdir(parents=True, exist_ok=True)
        num = sr.select_dtypes(include="number").columns
        sr[num] = sr[num].round(6)
        sr.to_csv(out_dir / f"{stem}.csv", index=False)
        head = [
            "# Top-N fastest / slowest race series by median WEATHER-CORRECTED v_j",
            "",
            "Estimand: per-series MEDIAN of the air-temp-corrected edition factor",
            "v_corr = v_j - b_w*(temp - temp_ref) -- every series evaluated at a",
            "COMMON reference temperature (net of typical weather; course + field",
            "stay in). 'Fastest if you catch a good-weather day', the inverse of",
            "q03's 'typical-weather-included' reading. Ranking is temp_ref-invariant;",
            "only the level depends on it. Median discards one freak edition (hence",
            "min-editions=3); selection by bootstrap rank stability P(top-N), medians",
            "recomputed inside each q07 replicate (carries v_j + b_w noise). Only",
            "air-temp-covered editions enter. min@3:00 = 180*(exp(median_v)-1).",
            "", f"slice: {slug}; min_n = {args.min_n}; min_editions = "
            f"{args.min_editions}; n_top = {args.n_top}; model = full_wxcorr", ""]
        (out_dir / f"{stem}.md").write_text("\n".join(head + md), encoding="utf-8")
        print(f"  wrote {out_dir / (stem + '.csv')}  ({len(sr)} rows)")
        print(f"  wrote {out_dir / (stem + '.md')}")
        n_done += 1

    if not n_done:
        raise SystemExit("no slices produced a result")

    print("\n".join(out_all))
    print()


if __name__ == "__main__":
    main()
