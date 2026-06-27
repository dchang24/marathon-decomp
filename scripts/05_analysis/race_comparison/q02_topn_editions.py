"""Top-N slowest / fastest race EDITIONS by v_j, with bootstrap rank stability.

ESTIMAND. v_j is the race factor under the beta=0 (bundling) gauge: the
expected log finish-time effect of that race edition for a reference athlete,
RELATIVE TO CONTEMPORANEOUS RACES. "Slowest" therefore means "slowest relative
to its era", which is the only cross-era reading the model supports. v_j
bundles course + that day's weather + field conditions -- a single hot edition
ranks on its own (see q03 for series-level rankings robust to one bad day).

SELECTION BY RANK STABILITY, NOT POINT RANKING. The point-estimate top-N is
contaminated by selection on noise (winner's curse): extreme v_j are extreme
partly because their sampling error drew extreme, and small-field races
dominate the tails. Instead, for every bootstrap replicate b (replicates share
the fit's gauge, so ranks are comparable across them):

    rank_slow_j^(b) = rank of v_j^(b), 1 = largest   (slowest)
    rank_fast_j^(b) = rank of v_j^(b), 1 = smallest  (fastest)
    P_slow_j = mean_b[ rank_slow_j^(b) <= N ]        (same for fast)

Reported sets per direction:
    headline  P >= p_headline (default 0.50)
    tie       p_tie <= P < p_headline (default 0.25) -- the "+1 or 2 more"
With no bootstrap for the slice, falls back to the point top-N (tier 'point').

Also per race: point rank, 95 pct CI of v (percentile across replicates), 95
pct rank interval, minutes-at-3:00 translation (180*(exp(v)-1) vs the average race).

The --min-n floor drops races with fewer than that many finishers IN THE SLICE
before any ranking (default 0 = keep all; rank-stability already demotes
small-noisy races, the floor is presentational). The floor, slice and model
are recorded in every output row.

Outputs, one dir per slice, filenames stamped with the CLI parameters so runs
with different settings never overwrite each other
(results/analysis/race_comparison/{slice}/):
    topn_editions_minn{min_n}_top{n_top}.csv   ALL eligible races of the slice
    topn_editions_minn{min_n}_top{n_top}.md    headline + tie tables
    + the same tables to stdout. The slice / min_n / n_top / model are also
    columns in the csv, so files from several slices concat cleanly.

Run::

    python scripts/05_analysis/race_comparison/q02_topn_editions.py
    python scripts/05_analysis/race_comparison/q02_topn_editions.py \
        --cohorts ALL Po10 --sexes M W B --min-n 25 --n-top 10
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
                  min_n: int, n_top: int, p_headline: float,
                  p_tie: float) -> pd.DataFrame | None:
    loaded = load_editions(slug, model, nutag, data_version)
    if loaded is None:
        return None
    df, fit_dir = loaded

    n_all = len(df)
    df = df[df["n_j"] >= min_n].reset_index(drop=True)
    n_excl = n_all - len(df)
    if n_excl:
        print(f"  {slug}: excluded {n_excl}/{n_all} races with n_j < {min_n}")
    if len(df) <= n_top:
        print(f"  [skip] {slug}: only {len(df)} eligible races (<= n_top)")
        return None

    v = df["v"].to_numpy(np.float64)
    B = load_boot_wide(fit_dir, df["race_id"].to_numpy(np.int64), v)
    has_boot = B is not None

    rank_slow, rank_fast = point_ranks(v)
    df["rank_slow"], df["rank_fast"] = rank_slow, rank_fast
    df["mult"] = np.exp(v)
    df["min_at_3h"] = v_to_min(v)

    p_slow = p_fast = None
    if has_boot:
        df["v_lo95"], df["v_hi95"] = pct(B, 2.5), pct(B, 97.5)
        rs, rf = rank_matrix(B), rank_matrix(-B)
        p_slow, p_fast = (rs <= n_top).mean(axis=0), (rf <= n_top).mean(axis=0)
        df["p_slow"], df["p_fast"] = p_slow, p_fast
        df["rank_slow_lo"], df["rank_slow_hi"] = pct(rs, 2.5), pct(rs, 97.5)
        df["rank_fast_lo"], df["rank_fast_hi"] = pct(rf, 2.5), pct(rf, 97.5)
        df["n_boot"] = B.shape[0]
    else:
        print(f"    [note] no bootstrap for {slug}; point-only tiers")
        for c in ("v_lo95", "v_hi95", "p_slow", "p_fast", "rank_slow_lo",
                  "rank_slow_hi", "rank_fast_lo", "rank_fast_hi"):
            df[c] = np.nan
        df["n_boot"] = 0

    df["tier_slow"] = tier_labels(p_slow, rank_slow, n_top, p_headline, p_tie)
    df["tier_fast"] = tier_labels(p_fast, rank_fast, n_top, p_headline, p_tie)

    df.insert(0, "slice", slug)
    df.insert(1, "model", f"{model}_{nutag}")
    df.insert(2, "min_n", min_n)
    df.insert(3, "n_top", n_top)
    # CSV row order: fastest -> slowest (most negative v first)
    df = df.sort_values("v").reset_index(drop=True)
    return df


# --------------------------------------------------------------------------- #
DIRS = {"slow": ("slowest", "rank_slow", "p_slow", "tier_slow",
                 "rank_slow_lo", "rank_slow_hi"),
        "fast": ("fastest", "rank_fast", "p_fast", "tier_fast",
                 "rank_fast_lo", "rank_fast_hi")}


def _fmt_rows(sub: pd.DataFrame, rank_c: str, p_c: str,
              rlo_c: str, rhi_c: str) -> list[list[str]]:
    rows = []
    for _, r in sub.iterrows():
        ci = (f"[{r['v_lo95']:+.4f}, {r['v_hi95']:+.4f}]"
              if pd.notna(r["v_lo95"]) else "-")
        rks = (f"[{r[rlo_c]:.0f}, {r[rhi_c]:.0f}]"
               if pd.notna(r[rlo_c]) else "-")
        rows.append([
            f"{int(r[rank_c])}", f"{r['series_key']} {int(r['year'])}",
            str(r["country"]), f"{int(r['n_j'])}", f"{r['v']:+.4f}", ci,
            f"{r['min_at_3h']:+.1f}",
            f"{r[p_c]:.2f}" if pd.notna(r[p_c]) else "-",
            rks, r["tier"] if "tier" in r else "",
        ])
    return rows


def report_slice(df: pd.DataFrame, n_top: int) -> tuple[list[str], list[str]]:
    """(markdown lines, stdout lines) for one slice's headline+tie tables."""
    md, out = [], []
    slug = df["slice"].iloc[0]
    hdr = ["rank", "race", "country", "n_j", "v_j", "v 95% CI",
           "min@3:00", "P(topN)", "rank 95%", "tier"]
    for key, (label, rank_c, p_c, tier_c, rlo_c, rhi_c) in DIRS.items():
        sub = df[df[tier_c] != ""].sort_values(rank_c).copy()
        sub["tier"] = sub[tier_c]
        title = f"{slug} -- {label} {n_top} (headline P>=0.5, tie P>=0.25)"
        md += [f"### {title}", "",
               "| " + " | ".join(hdr) + " |",
               "|" + "|".join(["---"] * len(hdr)) + "|"]
        out += ["", title,
                f"{'rank':>4} {'race':<28} {'ctry':>4} {'n_j':>6} {'v_j':>8} "
                f"{'min@3h':>7} {'P':>5} {'rank95':>10}  tier"]
        for row in _fmt_rows(sub, rank_c, p_c, rlo_c, rhi_c):
            md.append("| " + " | ".join(row) + " |")
            out.append(f"{row[0]:>4} {row[1]:<28} {row[2]:>4} {row[3]:>6} "
                       f"{row[4]:>8} {row[6]:>7} {row[7]:>5} {row[8]:>10}  {row[9]}")
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
                    help="drop races with fewer finishers in the slice (default 0 = keep all).")
    ap.add_argument("--p-headline", type=float, default=0.50)
    ap.add_argument("--p-tie", type=float, default=0.25)
    args = ap.parse_args()

    stem = f"topn_editions_minn{args.min_n}_top{args.n_top}"
    out_all, n_done = [], 0
    for cohort in args.cohorts:
        for sex in args.sexes:
            for mrc in args.mrcs:
                slug = f"{cohort}_{sex}_{args.window}_{mrc}"
                print(f"analysing {slug} ...")
                df = analyse_slice(slug, model=args.model, nutag=args.nutag,
                                   data_version=args.data_version,
                                   min_n=args.min_n, n_top=args.n_top,
                                   p_headline=args.p_headline, p_tie=args.p_tie)
                if df is None:
                    continue
                md, out = report_slice(df, args.n_top)
                out_all += out

                out_dir = OUT_ROOT / slug
                out_dir.mkdir(parents=True, exist_ok=True)
                num = df.select_dtypes(include="number").columns
                df[num] = df[num].round(6)
                df.to_csv(out_dir / f"{stem}.csv", index=False)
                head = [
                    "# Top-N slowest / fastest race editions by v_j",
                    "",
                    "Estimand: v_j under the beta=0 (bundling) gauge = race effect RELATIVE",
                    "TO CONTEMPORANEOUS RACES; it bundles course + that day's weather +",
                    "field conditions. Selection is by bootstrap rank stability P(top-N)",
                    "(headline P>=0.5, tie P>=0.25), which demotes small-field races whose",
                    "extreme point v_j is noise (winner's curse). min@3:00 =",
                    "180*(exp(v_j)-1) = minutes vs the average race for a 3:00:00 runner.",
                    "", f"slice: {slug}; min_n = {args.min_n}; n_top = {args.n_top}; "
                    f"model = {args.model}_{args.nutag}", ""]
                (out_dir / f"{stem}.md").write_text(
                    "\n".join(head + md), encoding="utf-8")
                print(f"  wrote {out_dir / (stem + '.csv')}  ({len(df)} rows)")
                print(f"  wrote {out_dir / (stem + '.md')}")
                n_done += 1

    if not n_done:
        raise SystemExit("no slices produced a result")

    print("\n".join(out_all))
    print()


if __name__ == "__main__":
    main()
