"""Top-N fastest / slowest race SERIES by the median edition v_j ("PB-chaser index").

ESTIMAND. Per series s with editions j in the slice, the index is

    m_s = median_j( v_j ),     over editions with n_j >= min_n,

with v_j under the beta=0 (bundling) gauge, i.e. era-relative by construction.
The median (not the mean) is the whole point: one freak edition (a hot Berlin
2025) is fully discarded once a series has >= 3 editions, with no hand-picked
exclusions -- the robust statistic formalizes "exclude the abnormal day"
without analyst discretion. The index INCLUDES each series' typical weather
(London's reliably cool April is legitimately part of why it is fast); a
net-of-weather "course index" variant is deferred to the covariate analysis.

m_s answers a runner's question directly: "relative to its contemporaries,
how fast does this race typically run?" (x180 = minutes at 3:00).

MIN EDITIONS (--min-editions, default 3). k = 3 is the smallest sample where
the median fully discards one outlying edition -- the property being bought.
Precision differences across series (k = 3 vs k = 11) are carried by the
bootstrap, not the floor. Unequal year coverage (2014-veterans vs post-COVID
arrivals) is handled by the estimand itself: every edition's v_j is already
relative to its own era. (Residual caveat: non-linear era effects -- supershoe
step, COVID -- are not removable by any gauge; check rank stability on a
common 2021-25 window as a sensitivity if needed.)

UNCERTAINTY AND TIES. For each bootstrap replicate b (replicates share the
fit's gauge), the series median is recomputed INSIDE the replicate -- this
preserves the within-replicate correlation between editions (shared athletes):

    m_s^(b) = median_j( v_j^(b) )
    rank_fast_s^(b) = rank of m_s^(b), 1 = smallest (fastest)
    P_fast_s = mean_b[ rank_fast_s^(b) <= N ]        (same for slow)

Reported per direction: headline set (P >= p_headline, default 0.50) and the
tie set (P >= p_tie, default 0.25) -- the "+1 or 2 more if there are ties".
No total order is claimed. With no bootstrap, falls back to the point top-N
(tier 'point').

The --min-n floor (finishers in the slice) is applied to EDITIONS before the
medians, so it changes which editions enter each series' index; floor, slice
and model are recorded in every output row.

Outputs, one dir per slice, filenames stamped with the CLI parameters so runs
with different settings never overwrite each other
(results/analysis/race_comparison/{slice}/):
    series_ranking_minn{min_n}_mined{min_editions}_top{n_top}.csv   all
                         eligible series of the slice
    series_ranking_minn{min_n}_mined{min_editions}_top{n_top}.md    headline
                         + tie tables
    + the same tables to stdout. The slice / min_n / min_editions / n_top /
    model are also columns in the csv, so files concat cleanly.

Run::

    python scripts/05_analysis/race_comparison/q03_series_ranking.py
    python scripts/05_analysis/race_comparison/q03_series_ranking.py \
        --cohorts ALL Po10 --sexes M W B --min-n 25 --min-editions 3 --n-top 10
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from race_common import (  # noqa: E402
    OUT_ROOT, v_to_min, load_boot_wide, load_editions,
    pct, point_ranks, rank_matrix, tier_labels,
)


# --------------------------------------------------------------------------- #
def analyse_slice(slug: str, *, model: str, nutag: str, data_version: str,
                  min_n: int, min_editions: int, n_top: int,
                  p_headline: float, p_tie: float) -> pd.DataFrame | None:
    loaded = load_editions(slug, model, nutag, data_version)
    if loaded is None:
        return None
    ed, fit_dir = loaded

    n_all = len(ed)
    ed = ed[ed["n_j"] >= min_n].reset_index(drop=True)
    if n_all - len(ed):
        print(f"  {slug}: excluded {n_all - len(ed)}/{n_all} editions with n_j < {min_n}")

    # series with enough editions after the floor
    sizes = ed.groupby("series_key").size()
    keep = sizes[sizes >= min_editions].index
    ed_k = ed[ed["series_key"].isin(keep)].reset_index(drop=True)
    print(f"  {slug}: {len(keep)} series with >= {min_editions} editions "
          f"({ed['series_key'].nunique() - len(keep)} dropped)")
    if len(keep) <= n_top:
        print(f"  [skip] {slug}: only {len(keep)} eligible series (<= n_top)")
        return None

    # point summary per series
    g = ed_k.groupby("series_key")
    sr = pd.DataFrame({
        "series_key": g.size().index,
        "country": g["country"].first().to_numpy(),
        "k_editions": g.size().to_numpy(int),
        "yr_min": g["year"].min().to_numpy(int),
        "yr_max": g["year"].max().to_numpy(int),
        "n_j_total": g["n_j"].sum().to_numpy(int),
        "median_v": g["v"].median().to_numpy(np.float64),
        "iqr_v": (g["v"].quantile(0.75) - g["v"].quantile(0.25)).to_numpy(np.float64),
    }).reset_index(drop=True)

    m = sr["median_v"].to_numpy(np.float64)
    sr["min_at_3h"] = v_to_min(m)
    rank_slow, rank_fast = point_ranks(m)
    sr["rank_slow"], sr["rank_fast"] = rank_slow, rank_fast

    # bootstrap: per-replicate medians, computed inside each replicate
    B = load_boot_wide(fit_dir, ed_k["race_id"].to_numpy(np.int64),
                       ed_k["v"].to_numpy(np.float64))
    has_boot = B is not None
    p_slow = p_fast = None
    if has_boot:
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
        print(f"    [note] no bootstrap for {slug}; point-only tiers")
        for c in ("m_lo95", "m_hi95", "p_slow", "p_fast", "rank_slow_lo",
                  "rank_slow_hi", "rank_fast_lo", "rank_fast_hi"):
            sr[c] = np.nan
        sr["n_boot"] = 0

    sr["tier_slow"] = tier_labels(p_slow, rank_slow, n_top, p_headline, p_tie)
    sr["tier_fast"] = tier_labels(p_fast, rank_fast, n_top, p_headline, p_tie)

    sr.insert(0, "slice", slug)
    sr.insert(1, "model", f"{model}_{nutag}")
    sr.insert(2, "min_n", min_n)
    sr.insert(3, "min_editions", min_editions)
    sr.insert(4, "n_top", n_top)
    # CSV row order: fastest -> slowest (most negative median_v first)
    sr = sr.sort_values("median_v").reset_index(drop=True)
    return sr


# --------------------------------------------------------------------------- #
DIRS = {"fast": ("fastest", "rank_fast", "p_fast", "tier_fast",
                 "rank_fast_lo", "rank_fast_hi"),
        "slow": ("slowest", "rank_slow", "p_slow", "tier_slow",
                 "rank_slow_lo", "rank_slow_hi")}


def report_slice(sr: pd.DataFrame, n_top: int) -> tuple[list[str], list[str]]:
    md, out = [], []
    slug = sr["slice"].iloc[0]
    hdr = ["rank", "series", "country", "k", "years", "median v",
           "m 95% CI", "min@3:00", "P(topN)", "rank 95%", "tier"]
    for key, (label, rank_c, p_c, tier_c, rlo_c, rhi_c) in DIRS.items():
        sub = sr[sr[tier_c] != ""].sort_values(rank_c)
        title = f"{slug} -- {label} {n_top} series (headline P>=0.5, tie P>=0.25)"
        md += [f"### {title}", "",
               "| " + " | ".join(hdr) + " |",
               "|" + "|".join(["---"] * len(hdr)) + "|"]
        out += ["", title,
                f"{'rank':>4} {'series':<26} {'ctry':>4} {'k':>3} {'years':>9} "
                f"{'med_v':>8} {'min@3h':>7} {'P':>5} {'rank95':>10}  tier"]
        for _, r in sub.iterrows():
            ci = (f"[{r['m_lo95']:+.4f}, {r['m_hi95']:+.4f}]"
                  if pd.notna(r["m_lo95"]) else "-")
            rks = (f"[{r[rlo_c]:.0f}, {r[rhi_c]:.0f}]"
                   if pd.notna(r[rlo_c]) else "-")
            yrs = f"{int(r['yr_min'])}-{int(r['yr_max'])%100:02d}"
            p_s = f"{r[p_c]:.2f}" if pd.notna(r[p_c]) else "-"
            row = [f"{int(r[rank_c])}", r["series_key"], str(r["country"]),
                   f"{int(r['k_editions'])}", yrs, f"{r['median_v']:+.4f}", ci,
                   f"{r['min_at_3h']:+.1f}", p_s, rks, r[tier_c]]
            md.append("| " + " | ".join(row) + " |")
            out.append(f"{row[0]:>4} {row[1]:<26} {row[2]:>4} {row[3]:>3} "
                       f"{row[4]:>9} {row[5]:>8} {row[7]:>7} {row[8]:>5} "
                       f"{row[9]:>10}  {row[10]}")
        md.append("")
    return md, out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohorts", nargs="+", default=["ALL", "Po10", "WA"])
    ap.add_argument("--sexes", nargs="+", default=["B", "M", "W"])
    ap.add_argument("--mrcs", nargs="+", default=["mrc2"])
    ap.add_argument("--window", default="14-25")
    ap.add_argument("--model", default="full", help="registry model tag (AxD = 'full').")
    ap.add_argument("--nutag", default="nu8p00")
    ap.add_argument("--data-version", default="race_results")
    ap.add_argument("--n-top", type=int, default=20)
    ap.add_argument("--min-n", type=int, default=0,
                    help="drop editions with fewer finishers in the slice before the medians.")
    ap.add_argument("--min-editions", type=int, default=3,
                    help="series needs >= this many (post-floor) editions; 3 = "
                         "smallest k where the median discards one outlier.")
    ap.add_argument("--p-headline", type=float, default=0.50)
    ap.add_argument("--p-tie", type=float, default=0.25)
    args = ap.parse_args()

    stem = (f"series_ranking_minn{args.min_n}"
            f"_mined{args.min_editions}_top{args.n_top}")
    out_all, n_done = [], 0
    for cohort in args.cohorts:
        for sex in args.sexes:
            for mrc in args.mrcs:
                slug = f"{cohort}_{sex}_{args.window}_{mrc}"
                print(f"analysing {slug} ...")
                sr = analyse_slice(slug, model=args.model, nutag=args.nutag,
                                   data_version=args.data_version,
                                   min_n=args.min_n,
                                   min_editions=args.min_editions,
                                   n_top=args.n_top,
                                   p_headline=args.p_headline, p_tie=args.p_tie)
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
                    "# Top-N fastest / slowest race series by median edition v_j",
                    "",
                    "Estimand ('PB-chaser index'): per-series MEDIAN of edition v_j under",
                    "the beta=0 (bundling) gauge -- era-relative by construction, and",
                    "INCLUDING the series' typical weather (a net-of-weather course index",
                    "is deferred to the covariate analysis). The median discards one freak",
                    "edition per series with no hand-picked exclusions (hence the",
                    "min-editions floor of 3). Selection is by bootstrap rank stability",
                    "P(top-N), medians recomputed inside each replicate; no total order is",
                    "claimed -- ties are the P>=0.25 set. min@3:00 = 180*(exp(median_v)-1) = minutes",
                    "vs the average race for a 3:00:00 runner.",
                    "", f"slice: {slug}; min_n = {args.min_n}; min_editions = "
                    f"{args.min_editions}; n_top = {args.n_top}; "
                    f"model = {args.model}_{args.nutag}", ""]
                (out_dir / f"{stem}.md").write_text(
                    "\n".join(head + md), encoding="utf-8")
                print(f"  wrote {out_dir / (stem + '.csv')}  ({len(sr)} rows)")
                print(f"  wrote {out_dir / (stem + '.md')}")
                n_done += 1

    if not n_done:
        raise SystemExit("no slices produced a result")

    print("\n".join(out_all))
    print()


if __name__ == "__main__":
    main()
