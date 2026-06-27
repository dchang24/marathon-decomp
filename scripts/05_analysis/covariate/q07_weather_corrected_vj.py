"""q07 - build a weather-corrected v_j (air temp only), point + bootstrap.

Removes the *day's air-temperature* contribution from each edition's race
factor, normalizing every edition to a reference temperature::

    v_corr_j = v_j  -  b_w * (temp_field_j - temp_ref)

where ``b_w`` is the **partial** air-temp slope from the joint, UNSTANDARDIZED
OLS ``v_j ~ temp_field + total_gain_m`` (course is a control only -- we subtract
weather alone, never course). ``temp_ref`` is the mean air temp over the
weather-covered editions of the slice; it only shifts the (gauge-arbitrary)
level, so it never changes any ranking. All ``v_j`` are the production full-model
factors under the G1 beta=0 APC gauge (via covariate_common); the correction is
applied in that gauge so it is consistent with the regression that produced
``b_w``. Higher v = slower = harder; v_corr = "how hard at reference weather,
still bundling course + field + the unexplained part."

JOINT BOOTSTRAP. ``b_w`` is itself noisy. For every fit bootstrap replicate b
(athlete resampling; the covariates are fixed race attributes) we re-gauge
v^(b), refit the joint OLS on v^(b) to get b_w^(b), and correct:
``v_corr_j^(b) = v^(b)_j - b_w^(b)*(temp_j - temp_ref)``. The downstream ranker
(race_comparison/q04) selects on these replicates, so its rank-stability
P(top-N) propagates BOTH v_j sampling noise and b_w uncertainty.

SAMPLES. ``b_w`` is estimated on the joint complete-case set (v, temp, gain all
present). The correction -- and every output row -- covers the editions with
air temp present (``covered``); a temp-present edition that lacks elevation is
still corrected (correction needs temp only) but did not enter the b_w fit.

Outputs (results/analysis/covariate/07_weather_corrected_vj/):
  weather_corrected_vj__6slices.parquet       one row per (slice, covered race):
      race_id, n_j, v (gauged point), temp_field, b_w (slice scalar), temp_ref,
      v_corr (point). Plus a per-slice b_w bootstrap SD column for reference.
  weather_corrected_vj_boot__6slices.parquet  one row per (slice, run_id, race):
      v_corr replicate (only when a bootstrap exists for the slice).
  + a per-slice summary to stdout (b_w point [boot sd], temp_ref, coverage).

Run after q01 (needs covariate_merged.parquet) ::

    python scripts/05_analysis/covariate/q07_weather_corrected_vj.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import covariate_common as C

SUBDIR = "07_weather_corrected_vj"
WEATHER = "temp_field"          # air temperature (C); q02 selection winner
COURSE = C.COURSE_PREDICTOR     # total_gain_m -- control only, never subtracted


def _ols_beta(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Plain OLS coefficient vector (no intercept handling; X has the column)."""
    return np.linalg.solve(X.T @ X, X.T @ y)


def correct_slice(merged: pd.DataFrame, slice_name: str
                  ) -> tuple[pd.DataFrame, pd.DataFrame | None, dict]:
    """Point + (optional) bootstrap weather-corrected v_j for one slice.

    Returns (point_df, boot_df_or_None, summary_dict). All covariates are taken
    from the merged table (keyed by global race_id); the point v_j is the gauged
    `v_{slice}` column; the bootstrap matrix comes re-gauged from
    covariate_common and is aligned to the merged covariates by race_id.
    """
    vcol, ncol = f"v_{slice_name}", f"n_{slice_name}"
    midx = merged.set_index("race_id")

    # ---- point: regress on joint complete cases, correct the temp-covered set
    in_slice = merged[vcol].notna()
    d = merged.loc[in_slice, ["race_id", vcol, ncol, WEATHER, COURSE]].copy()
    v = d[vcol].to_numpy(float)
    temp = d[WEATHER].to_numpy(float)
    gain = d[COURSE].to_numpy(float)
    cc = np.isfinite(v) & np.isfinite(temp) & np.isfinite(gain)   # b_w fit set
    cov = np.isfinite(temp)                                       # correctable
    if cc.sum() < 20 or cov.sum() == 0:
        raise SystemExit(f"{slice_name}: too few complete cases (cc={int(cc.sum())})")

    temp_ref = float(np.mean(temp[cov]))
    Xcc = np.column_stack([np.ones(int(cc.sum())), temp[cc], gain[cc]])
    b_w = float(_ols_beta(v[cc], Xcc)[1])
    v_corr = v - b_w * (temp - temp_ref)

    point = pd.DataFrame({
        "slice": slice_name,
        "race_id": d["race_id"].to_numpy()[cov],
        "n_j": d[ncol].to_numpy()[cov].astype(int),
        "v": v[cov],
        "temp_field": temp[cov],
        "b_w": b_w,
        "temp_ref": temp_ref,
        "v_corr": v_corr[cov],
    })

    summary = dict(slice=slice_name, n_slice=int(in_slice.sum()),
                   n_cc=int(cc.sum()), n_cov=int(cov.sum()),
                   b_w=b_w, temp_ref=temp_ref, b_w_boot_sd=np.nan, R=0)

    # ---- bootstrap: refit b_w per replicate, correct, emit long rows
    boot = C.load_vj_boot_wide(slice_name)
    boot_df = None
    if boot is not None:
        _fdir, race_ids, _t, Bg = boot          # Bg: (R, J) re-gauged replicates
        tempB = midx.reindex(race_ids)[WEATHER].to_numpy(float)
        gainB = midx.reindex(race_ids)[COURSE].to_numpy(float)
        finite = np.isfinite(Bg).all(axis=0)    # races present in every replicate
        ccB = np.isfinite(tempB) & np.isfinite(gainB) & finite
        covB = np.isfinite(tempB) & finite
        Xb = np.column_stack([np.ones(int(ccB.sum())), tempB[ccB], gainB[ccB]])
        XtXinv = np.linalg.inv(Xb.T @ Xb)
        coef = XtXinv @ (Xb.T @ Bg[:, ccB].T)   # (3, R)
        bw_rep = coef[1]                          # (R,) per-replicate temp slope
        # v_corr^(b)_j = v^(b)_j - bw^(b) * (temp_j - temp_ref), covered races only
        Vc = Bg[:, covB] - np.outer(bw_rep, tempB[covB] - temp_ref)
        R = Bg.shape[0]
        rid_cov = race_ids[covB]
        boot_df = pd.DataFrame({
            "slice": slice_name,
            "run_id": np.repeat(np.arange(1, R + 1), rid_cov.size),
            "race_id": np.tile(rid_cov, R),
            "v_corr": Vc.reshape(-1),
        })
        summary["b_w_boot_sd"] = float(np.std(bw_rep, ddof=1))
        summary["R"] = R
        # carry the bootstrap b_w SD onto the point table for reference
        point["b_w_boot_sd"] = summary["b_w_boot_sd"]
    else:
        point["b_w_boot_sd"] = np.nan

    return point, boot_df, summary


def main() -> None:
    merged = pd.read_parquet(C.MERGED_PATH)

    pts, boots, summ = [], [], []
    for s in C.SLICE_ORDER:
        point, boot_df, summary = correct_slice(merged, s)
        pts.append(point)
        if boot_df is not None:
            boots.append(boot_df)
        summ.append(summary)

    tag = C.ALL_SLICES_TAG
    point_path = C.out_path(SUBDIR, "weather_corrected_vj", tag, "parquet")
    boot_path = C.out_path(SUBDIR, "weather_corrected_vj_boot", tag, "parquet")
    point_all = pd.concat(pts, ignore_index=True)
    point_all.to_parquet(point_path, index=False)
    if boots:
        pd.concat(boots, ignore_index=True).to_parquet(boot_path, index=False)

    print("=" * 86)
    print("q07  weather-corrected v_j   v_corr = v - b_w*(temp - temp_ref)   (air temp only)")
    print("     b_w = partial temp slope from joint OLS v ~ temp + total_gain (UNSTANDARDIZED)")
    print("=" * 86)
    print(f"  {'slice':8s} {'n_slice':>7s} {'n_cc':>5s} {'n_cov':>5s} "
          f"{'b_w (log/C)':>14s} {'[boot sd]':>11s} {'temp_ref':>9s}")
    for r in summ:
        bsd = f"{r['b_w_boot_sd']:.5f}" if np.isfinite(r["b_w_boot_sd"]) else "-"
        # log/C -> min @ 3:00 per +1 C, a tangible sense of the slope size
        min_per_c = 180.0 * r["b_w"]
        print(f"  {r['slice']:8s} {r['n_slice']:7d} {r['n_cc']:5d} {r['n_cov']:5d} "
              f"{r['b_w']:+14.5f} {bsd:>11s} {r['temp_ref']:9.2f}"
              f"   ({min_per_c:+.2f} min@3:00 per +1C)")

    print(f"\n[write] {point_path}")
    if boots:
        print(f"[write] {boot_path}")
    else:
        print("[note] no bootstrap found for any slice; point-only correction written")


if __name__ == "__main__":
    main()
