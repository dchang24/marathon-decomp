"""q04 - does the v_j ~ weather + course regression depend on outliers or on
       v_j athlete-sampling uncertainty?

For each of the six slices we fit the same joint standardized model
    v_j ~ z(temp_field) + z(total_gain_m)
under a 2x2 design (four fits per slice):

    OLS         unweighted, Gaussian          (= q03 point fit)
    WLS         weight 1/var_j, Gaussian      (var_j = bootstrap SD^2 of v_j)
    Robust      unweighted, Huber c=1.345     (down-weight outlier races)
    Robust+WLS  weight 1/var_j, Huber

`var_j` is the across-bootstrap-replicate variance of the beta=0-gauged v_j
(athlete resampling) -> WLS trusts well-determined races more. Comparing the
four standardized betas shows how much the conclusion leans on a few outlier
races (OLS vs Robust) and on the precision of v_j (unweighted vs weighted).

The simple metrics carry no v_j bootstrap uncertainty, so for them we report only
the outlier axis (OLS vs Robust, unweighted) as a secondary check.

Console + results/analysis/covariate/04_regression_robustness/
regression_robustness__6slices.{parquet,csv}. The flagged outlier races (lowest
Huber weight) print per slice and save to regression_outliers__6slices.csv.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import covariate_common as C

SUBDIR = "04_regression_robustness"
WEATHER = "temp_field"      # q02 winner; the headline weather predictor
COURSE = C.COURSE_PREDICTOR
FITS = ["OLS", "WLS", "Robust", "Robust+WLS"]


def transform(df: pd.DataFrame, col: str) -> np.ndarray:
    x = df[col].to_numpy(float)
    if col.endswith("_sec"):
        return np.where(x > 0, np.log(x), np.nan)
    return x


def four_fits(y: np.ndarray, Xz: np.ndarray, var: np.ndarray | None) -> dict:
    """Return {fit_kind: (beta, se, r2, W)} for the available fits.

    `var` is the per-row v_j sampling variance (None -> no weighted fits).
    Xz is the standardized predictor block (no intercept); y is standardized.
    """
    n = len(y)
    X = np.column_stack([np.ones(n), Xz])
    out = {}
    pw = None
    if var is not None:
        pw = 1.0 / np.maximum(var, 1e-12)
        pw = pw / pw.mean()                  # normalize so weighted R2 ~ comparable
    out["OLS"] = C.robust_wls(X, y, None, huber=False)
    out["Robust"] = C.robust_wls(X, y, None, huber=True)
    if pw is not None:
        out["WLS"] = C.robust_wls(X, y, pw, huber=False)
        out["Robust+WLS"] = C.robust_wls(X, y, pw, huber=True)
    return out


def main() -> None:
    df = pd.read_parquet(C.MERGED_PATH)

    # attach per-race v_j sampling SD for each slice
    for s in C.SLICE_ORDER:
        sd = C.load_vj_boot_sd(s)
        if sd is not None:
            df = df.merge(sd, on="race_id", how="left")

    rows, outlier_rows = [], []

    # ---------------- (A) the six slices: full 2x2 ---------------- #
    print("=" * 92)
    print("A. v_j ~ temp + total_gain   4 fits per slice  (standardized betas; z in [])")
    print("   OLS=q03 point | WLS=1/var_j | Robust=Huber | Robust+WLS")
    print("=" * 92)
    for s in C.SLICE_ORDER:
        ycol, vsdcol = f"v_{s}", f"vsd_{s}"
        y_raw = df[ycol].to_numpy(float)
        var = df[vsdcol].to_numpy(float) ** 2 if vsdcol in df else None
        m = np.isfinite(y_raw) & np.isfinite(df[WEATHER].to_numpy(float)) \
            & np.isfinite(df[COURSE].to_numpy(float))
        if var is not None:
            m &= np.isfinite(df[vsdcol].to_numpy(float))
        n = int(m.sum())
        yz = C.standardize(y_raw[m])
        Xz = np.column_stack([C.standardize(df[WEATHER].to_numpy(float)[m]),
                              C.standardize(df[COURSE].to_numpy(float)[m])])
        fits = four_fits(yz, Xz, var[m] if var is not None else None)

        print(f"\n  {s}  (n={n})")
        print(f"    {'fit':12s} {'b_temp[z]':>14s} {'b_gain[z]':>14s} {'R2':>7s} "
              f"{'n_down':>7s}")
        for fk in FITS:
            if fk not in fits:
                continue
            beta, se, r2, W = fits[fk]
            z = beta / se
            # Huber-downweight count is only interpretable for the pure-Huber fit;
            # under WLS the weight spread is dominated by the 1/var prior.
            n_down = int((W < 0.5).sum()) if fk == "Robust" else -1
            nd = f"{n_down:d}" if n_down >= 0 else "."
            print(f"    {fk:12s} {beta[1]:+.3f}[{z[1]:+5.1f}] {beta[2]:+.3f}[{z[2]:+5.1f}]"
                  f" {r2:7.3f} {nd:>7s}")
            rows.append(dict(target=ycol, kind="slice", fit=fk, n=n,
                             b_temp=beta[1], z_temp=z[1], b_gain=beta[2],
                             z_gain=z[2], r2=r2))
        # outliers from the Robust (unweighted) fit
        beta_r, _, _, W_r = fits["Robust"]
        rid = df["race_id"].to_numpy()[m]
        ser = df["series_key"].to_numpy()[m]
        yr = df["year"].to_numpy()[m]
        resid = yz - np.column_stack([np.ones(n), Xz]) @ beta_r
        order = np.argsort(W_r)            # smallest weight = most down-weighted
        flagged = order[W_r[order] < 0.7][:6]
        if len(flagged):
            print(f"    outliers (Huber w<0.7): "
                  + ", ".join(f"{ser[i]}{int(yr[i])}(w={W_r[i]:.2f})" for i in flagged))
        for i in flagged:
            outlier_rows.append(dict(slice=s, race_id=int(rid[i]), series_key=ser[i],
                                     year=int(yr[i]), huber_w=float(W_r[i]),
                                     std_resid=float(resid[i])))

    # ---- per-slice deltas: how much do betas move? ----
    print("\n" + "-" * 92)
    print("  beta_temp shift vs OLS (how much each axis matters)")
    print(f"    {'slice':10s} {'OLS':>8s} {'WLS-OLS':>9s} {'Rob-OLS':>9s} {'RobW-OLS':>9s}")
    rr = pd.DataFrame(rows)
    for s in C.SLICE_ORDER:
        sub = rr[(rr.target == f"v_{s}")].set_index("fit")
        if "OLS" not in sub.index:
            continue
        b0 = sub.loc["OLS", "b_temp"]
        d = {fk: (sub.loc[fk, "b_temp"] - b0 if fk in sub.index else np.nan) for fk in FITS}
        print(f"    {s:10s} {b0:+8.3f} {d['WLS']:+9.3f} {d['Robust']:+9.3f} "
              f"{d['Robust+WLS']:+9.3f}")

    # ---------------- (B) simple metrics: outlier axis only ---------------- #
    print("\n" + "=" * 92)
    print("B. Simple metrics ~ temp + total_gain   OLS vs Robust (outlier sensitivity)")
    print("   (no v_j uncertainty -> unweighted only; |b_temp shift| ranked, top 12)")
    print("=" * 92)
    metric_cols = [c for c in df.columns if c.endswith("_sec") or c.endswith("_pct")]
    mrows = []
    for col in metric_cols:
        y_raw = transform(df, col)
        m = np.isfinite(y_raw) & np.isfinite(df[WEATHER].to_numpy(float)) \
            & np.isfinite(df[COURSE].to_numpy(float))
        n = int(m.sum())
        if n < 20:
            continue
        yz = C.standardize(y_raw[m])
        Xz = np.column_stack([C.standardize(df[WEATHER].to_numpy(float)[m]),
                              C.standardize(df[COURSE].to_numpy(float)[m])])
        X = np.column_stack([np.ones(n), Xz])
        bo, _, r2o, _ = C.robust_wls(X, yz, None, huber=False)
        br, _, r2r, _ = C.robust_wls(X, yz, None, huber=True)
        mrows.append(dict(target=col, n=n, b_temp_ols=bo[1], b_temp_rob=br[1],
                          shift=br[1] - bo[1], r2_ols=r2o, r2_rob=r2r))
        for fk, (b, r2) in [("OLS", (bo, r2o)), ("Robust", (br, r2r))]:
            rows.append(dict(target=col, kind="metric", fit=fk, n=n,
                             b_temp=b[1], z_temp=np.nan, b_gain=b[2],
                             z_gain=np.nan, r2=r2))
    md = pd.DataFrame(mrows)
    md["abs_shift"] = md["shift"].abs()
    top = md.sort_values("abs_shift", ascending=False).head(12)
    with pd.option_context("display.float_format", lambda v: f"{v:+.3f}"):
        print(top[["target", "n", "b_temp_ols", "b_temp_rob", "shift",
                   "r2_ols", "r2_rob"]].to_string(index=False))

    out = pd.DataFrame(rows)
    tag = C.ALL_SLICES_TAG
    pq = C.out_path(SUBDIR, "regression_robustness", tag, "parquet")
    out.to_parquet(pq, index=False)
    out.to_csv(C.out_path(SUBDIR, "regression_robustness", tag, "csv"), index=False)
    outliers_csv = C.out_path(SUBDIR, "regression_outliers", tag, "csv")
    pd.DataFrame(outlier_rows).to_csv(outliers_csv, index=False)
    print(f"\n[write] {pq}")
    print(f"[write] {outliers_csv}")


if __name__ == "__main__":
    main()
