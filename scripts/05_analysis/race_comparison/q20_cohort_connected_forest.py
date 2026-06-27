"""q20 - compute + EXPORT the processed data behind the cohort connected forest
figure (p20 reads these; q20 does the numbers, p20 only draws).

For each slice (ALL_B, Po10_B) it assembles one tidy row per series with every
number the figure shows:

  series info   country, k_editions, yr_min, yr_max, n_men, n_women, n_total
  ranking       rank (1 = fastest by the BOOTSTRAP MEDIAN), tier (dominance,
                thr 0.9, recomputed in bootstrap-median order so bands match
                the plotted order)
  point         point_median_v          (q05 point-fit series median)
  bootstrap     boot_min, boot_p025, boot_p25, boot_p50 (= ranking metric),
                boot_p75, boot_p975, boot_max, boot_whislo, boot_whishi
                (the exact boxplot stats: box = p25/p50/p75, matplotlib 1.5*IQR
                whiskers = whislo/whishi); +boot_median_min3h = 180*(exp(boot_p50)-1)
  sex (colour)  sex_v_men, sex_v_women      per-sex series medians (q04 MvsW,
                women re-pinned into men's gauge), sex_contrast_r /
                sex_contrast_min3h = the gauge-free within-cohort M-vs-W
                contrast, sex_flag

Plus a connector table (one row per shared series) with the across-cohort
ALL-vs-Po10 gauge-free contrast + flag that drives the thick connectors.

Ranking metric is now the BOOTSTRAP MEDIAN (boot_p50), so the boxplot centre
line is monotone with row order (the old point-vs-bootstrap mismatch is gone).

Inputs (run q05 ALL_B+Po10_B and q04 ALLvsPo10_B / MvsW_ALL / MvsW_Po10 first):
    {ALL_B,Po10_B}/series_dominance.csv, series_boot_medians.parquet
    cross_slice/ALLvsPo10_B/cross_movers_series.csv
    cross_slice/MvsW_{ALL,Po10}/cross_movers_series.csv
Outputs (results/analysis/race_comparison/cohort_forest/):
    cohort_forest_ALL_B.csv, cohort_forest_Po10_B.csv      per-series numbers
    cohort_forest_ALL_B_boot.parquet, ..._Po10_B_boot.parquet  raw bootstrap
                                                           medians (for p20 markers)
    cohort_forest_connectors.csv                           shared-series connectors

Run::

    python scripts/05_analysis/race_comparison/q20_cohort_connected_forest.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from race_common import OUT_ROOT, v_to_min  # noqa: E402

CF_DIR = OUT_ROOT / "cohort_forest"
SLICES = {"ALL": "ALL_B_14-25_mrc2", "Po10": "Po10_B_14-25_mrc2"}
SEX_PAIR = {"ALL": "MvsW_ALL_14-25_mrc2", "Po10": "MvsW_Po10_14-25_mrc2"}
COHORT_PAIR = "ALLvsPo10_B_14-25_mrc2"


def box_stats(data: np.ndarray) -> dict:
    p = np.percentile(data, [2.5, 25, 50, 75, 97.5])
    q1, q3 = p[1], p[3]
    iqr = q3 - q1
    hi, lo = q3 + 1.5 * iqr, q1 - 1.5 * iqr
    whishi = float(data[data <= hi].max())
    whislo = float(data[data >= lo].min())
    return dict(boot_min=float(data.min()), boot_p025=p[0], boot_p25=q1,
                boot_p50=p[2], boot_p75=q3, boot_p975=p[4],
                boot_max=float(data.max()), boot_whislo=whislo, boot_whishi=whishi)


def dominance_tiers(M: np.ndarray, order: np.ndarray, thr: float) -> np.ndarray:
    """Leader-sweep tiers in `order` (fastest-first); 1 = fastest tier."""
    S = M.shape[1]
    D = np.empty((S, S))
    for s in range(S):
        D[s] = (M[:, [s]] < M).mean(axis=0)
    np.fill_diagonal(D, np.nan)
    tiers = np.empty(S, int)
    t, leader = 0, order[0]
    for k, s in enumerate(order):
        if k > 0 and D[leader, s] >= thr:
            t += 1
            leader = s
        tiers[s] = t
    return tiers + 1


def load_sex(cohort: str) -> pd.DataFrame:
    p = OUT_ROOT / "cross_slice" / SEX_PAIR[cohort] / "cross_movers_series.csv"
    if not p.is_file():
        print(f"  [warn] missing {p}; no sex contrast for {cohort}")
        return pd.DataFrame(columns=["series_key", "sex_v_men", "sex_v_women",
                                     "sex_contrast_r", "sex_contrast_min3h", "sex_flag"])
    d = pd.read_csv(p)
    return pd.DataFrame({
        "series_key": d["series_key"],
        "sex_v_men": d["median_vA"], "sex_v_women": d["median_vBstar"],
        "sex_contrast_r": d["contrast_r"], "sex_contrast_min3h": d["r_min_at_3h"],
        "sex_flag": d["flag"].astype(bool),
    })


def build_slice(cohort: str, thr: float):
    slug = SLICES[cohort]
    rank = pd.read_csv(OUT_ROOT / slug / "series_dominance.csv")
    boot = pd.read_parquet(OUT_ROOT / slug / "series_boot_medians.parquet")
    M = np.column_stack([boot[s].to_numpy(np.float64) for s in rank["series_key"]])

    stats = pd.DataFrame([box_stats(M[:, i]) for i in range(M.shape[1])])
    df = pd.concat([rank[["series_key", "country", "k_editions", "yr_min",
                          "yr_max", "n_men", "n_women", "median_v"]]
                    .reset_index(drop=True), stats], axis=1)
    df = df.rename(columns={"median_v": "point_median_v"})
    df["n_total"] = df["n_men"] + df["n_women"]
    df["boot_median_min3h"] = v_to_min(df["boot_p50"])

    order = np.argsort(df["boot_p50"].to_numpy())          # fastest first
    df["tier"] = dominance_tiers(M, order, thr)

    df = df.merge(load_sex(cohort), on="series_key", how="left")
    df = df.sort_values("boot_p50").reset_index(drop=True)
    df.insert(0, "rank", np.arange(1, len(df) + 1))
    df.insert(0, "slice", slug)

    cols = ["slice", "rank", "series_key", "country", "k_editions", "yr_min",
            "yr_max", "tier", "n_men", "n_women", "n_total", "point_median_v",
            "boot_p50", "boot_median_min3h", "boot_min", "boot_p025", "boot_p25",
            "boot_p75", "boot_p975", "boot_max", "boot_whislo", "boot_whishi",
            "sex_v_men", "sex_v_women", "sex_contrast_r", "sex_contrast_min3h",
            "sex_flag"]
    return df[cols], boot, slug


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--thr", type=float, default=0.9,
                    help="dominance threshold for the tier sweep.")
    args = ap.parse_args()

    CF_DIR.mkdir(parents=True, exist_ok=True)
    for cohort in ("ALL", "Po10"):
        df, boot, slug = build_slice(cohort, args.thr)
        num = df.select_dtypes(include="number").columns
        df[num] = df[num].round(6)
        df.to_csv(CF_DIR / f"cohort_forest_{cohort}_B.csv", index=False)
        boot.to_parquet(CF_DIR / f"cohort_forest_{cohort}_B_boot.parquet", index=False)
        print(f"  wrote cohort_forest_{cohort}_B.csv  ({len(df)} series, "
              f"{df['tier'].max()} tiers) + _boot.parquet")

    # connectors: across-cohort ALL-vs-Po10 contrast + flag (shared series)
    cp = OUT_ROOT / "cross_slice" / COHORT_PAIR / "cross_movers_series.csv"
    if cp.is_file():
        c = pd.read_csv(cp)
        conn = pd.DataFrame({
            "series_key": c["series_key"], "country": c["country"],
            "k_editions": c["k_editions"],
            "cohort_contrast_r": c["contrast_r"],
            "cohort_contrast_min3h": c["r_min_at_3h"],
            "cohort_z": c["z"], "cohort_flag": c["flag"].astype(bool),
            "rank_ALL": c["rank_A"], "rank_Po10": c["rank_B"], "d_rank": c["d_rank"],
        }).sort_values("cohort_z", ascending=False, na_position="last")
        num = conn.select_dtypes(include="number").columns
        conn[num] = conn[num].round(6)
        conn.to_csv(CF_DIR / "cohort_forest_connectors.csv", index=False)
        print(f"  wrote cohort_forest_connectors.csv  ({len(conn)} shared, "
              f"{int(conn['cohort_flag'].sum())} flagged)")
    else:
        print(f"  [warn] missing {cp}; no connector table")


if __name__ == "__main__":
    main()
