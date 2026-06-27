"""q01 - build the merged per-race covariate table + coverage audit.

Pulls the production **full**-model v_j for the six slices (ALL/Po10 x M/W/B),
beta=0-gauges each, and joins:
  * weather    (per competition_id): all WBGT/temp/wind/precip features
  * elevation  (per series_key, year): net / total gain (+ |net|)
  * simple metrics (per competition_id): every *_sec / *_pct in
    data/misc/race_stats_summary.csv

Outputs (results/analysis/covariate/01_merge/):
  * covariate_merged__6slices.parquet  - full merged table (consumed by q02-q07 + p01)
  * series_year_covariates.csv         - focused per series+year lookup of air temp,
    WBGT (field/max) and course elevation gain (human reference; also read by
    scripts/visualizations/race_dashboard.py for hover info)
plus a coverage audit to the console. Build + audit only.
"""
from __future__ import annotations

import covariate_common as C


def main() -> None:
    df, metric_cols = C.build_merged()

    # coverage is judged on the ALL_B race set (the superset = every race)
    base = df[df[f"v_{C.VAR_SELECT_SLICE}"].notna()]
    n = len(base)
    has_temp = base["temp_field"].notna()
    has_wbgt = base["wbgt_max"].notna()
    has_el = base["total_gain_m"].notna()
    print(f"\n[coverage] {n} races in {C.VAR_SELECT_SLICE} (superset slice)")
    print(f"  air temp      : {has_temp.sum():4d} ({100 * has_temp.mean():.0f}%)")
    print(f"  WBGT max      : {has_wbgt.sum():4d} ({100 * has_wbgt.mean():.0f}%)")
    print(f"  total gain    : {has_el.sum():4d} ({100 * has_el.mean():.0f}%)")
    print(f"  weather+course: {(has_temp & has_el).sum():4d} "
          f"({100 * (has_temp & has_el).mean():.0f}%)")
    print(f"  simple metrics: {len(metric_cols)} columns "
          f"({', '.join(metric_cols[:6])}, ...)")

    print("\n[per-slice v_j coverage]")
    for s in C.SLICE_ORDER:
        print(f"  {s:8s}: {df[f'v_{s}'].notna().sum():4d} races")

    out = C.MERGED_PATH
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out, index=False)
    print(f"\n[write] {out}  ({len(df.columns)} cols, {len(df)} rows)")

    # focused per series+year covariate lookup (air temp + WBGT + elevation gain)
    cov = (df[C.SERIES_YEAR_COV_COLS]
           .sort_values(["series_key", "year"])
           .reset_index(drop=True))
    cov.to_csv(C.SERIES_YEAR_COV_PATH, index=False)
    n_temp = cov["temp_field"].notna().sum()
    n_gain = cov["total_gain_m"].notna().sum()
    print(f"[write] {C.SERIES_YEAR_COV_PATH}  ({len(cov)} races; "
          f"temp {n_temp}, gain {n_gain})")

    missing = base[~has_temp & ~has_el][["race_id", "series_key", "year"]]
    if len(missing):
        print(f"\n[no-covariate races] {len(missing)}:")
        print(missing.to_string(index=False, max_rows=40))


if __name__ == "__main__":
    main()
