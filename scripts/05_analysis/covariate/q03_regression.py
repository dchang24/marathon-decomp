"""q03 - regress each target on weather + course jointly, across slices.

Two-predictor standardized OLS:  target ~ weather + total_gain_m, run once per
weather choice (the q02 winner `temp_field`, plus the principled `wbgt_max`):
  * temp model : y ~ z(temp_field) + z(total_gain_m)
  * WBGT model : y ~ z(wbgt_max)   + z(total_gain_m)

Targets (two families):
  (A) full-model v_j for each of the six slices (ALL/Po10 x M/W/B)
  (B) every naive 'simple' metric in race_stats_summary.csv (field-contaminated
      proxies: finish-time percentiles, top-N elites, sub-N / BQ shares)

For each target x weather-set we report the joint R2 (bootstrap 95% CI), and the
standardized partial beta + z of each predictor. Finish-time *_sec metrics are
log-transformed (multiplicative, like v_j); *_pct rise when a race is EASIER.

Headline: how much of each slice's v_j is explained by objective conditions, and
how that compares against the simple metrics. Console +
results/analysis/covariate/03_regression/regression__6slices.{parquet,csv}.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import covariate_common as C

SUBDIR = "03_regression"
COURSE = C.COURSE_PREDICTOR


def transform(df: pd.DataFrame, col: str) -> np.ndarray:
    """log for finish-time seconds (multiplicative, like v_j); identity else."""
    x = df[col].to_numpy(float)
    if col.endswith("_sec"):
        return np.where(x > 0, np.log(x), np.nan)
    return x


def _r2(X, y):
    b = np.linalg.solve(X.T @ X, X.T @ y)
    resid = y - X @ b
    return 1 - (resid @ resid) / ((y - y.mean()) ** 2).sum()


def boot_r2_ci(X, y, *, B: int = 2000, seed: int = 0) -> tuple[float, float, float]:
    """Nonparametric case-resampling CI + SD for R2. Returns (lo, hi, sd)."""
    rng = np.random.default_rng(seed)
    n = len(y)
    draws = np.empty(B)
    for b in range(B):
        idx = rng.integers(0, n, n)
        draws[b] = _r2(X[idx], y[idx])
    return (float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5)),
            float(draws.std(ddof=1)))


def fit_one(df: pd.DataFrame, ycol: str, wcol: str, *, boot: int = 2000) -> dict | None:
    y = transform(df, ycol)
    preds = [wcol, COURSE]
    m = np.isfinite(y)
    for p in preds:
        m &= np.isfinite(df[p].to_numpy(float))
    n = int(m.sum())
    if n < 20:
        return None
    yz = C.standardize(y[m])
    Xz = np.column_stack([C.standardize(df[p].to_numpy(float)[m]) for p in preds])
    X = np.column_stack([np.ones(n), Xz])
    beta, se, r2 = C.ols(X, yz)
    z = beta / se
    zc = 1.959963985  # 95% normal
    r2_lo, r2_hi, r2_sd = boot_r2_ci(X, yz, B=boot)
    return dict(n=n, r2=r2, r2_sd=r2_sd, r2_lo=r2_lo, r2_hi=r2_hi,
                beta_weather=beta[1], se_weather=se[1], z_weather=z[1],
                beta_weather_lo=beta[1] - zc * se[1], beta_weather_hi=beta[1] + zc * se[1],
                beta_course=beta[2], se_course=se[2], z_course=z[2],
                beta_course_lo=beta[2] - zc * se[2], beta_course_hi=beta[2] + zc * se[2])


def main() -> None:
    df = pd.read_parquet(C.MERGED_PATH)
    targets = [(f"v_{s}", "slice") for s in C.SLICE_ORDER]
    metric_cols = [c for c in df.columns if c.endswith("_sec") or c.endswith("_pct")]
    targets += [(c, "metric") for c in metric_cols]

    rows = []
    for ycol, kind in targets:
        for wcol, wlab in C.WEATHER_SETS:
            res = fit_one(df, ycol, wcol)
            if res is None:
                continue
            rows.append(dict(target=ycol, kind=kind, weather_set=wlab, **res))
    out = pd.DataFrame(rows)
    pq = C.out_path(SUBDIR, "regression", C.ALL_SLICES_TAG, "parquet")
    out.to_parquet(pq, index=False)
    out.to_csv(C.out_path(SUBDIR, "regression", C.ALL_SLICES_TAG, "csv"), index=False)

    # ---- (A) the six slices head to head ----
    print("=" * 86)
    print("A. v_j ~ weather + total_gain  (joint 2-predictor standardized OLS)")
    print("   how much of each slice's v_j is explained by objective conditions?")
    print("=" * 86)
    sl = out[out.kind == "slice"]
    for wlab in [w[1] for w in C.WEATHER_SETS]:
        sub = sl[sl.weather_set == wlab].set_index("target")
        print(f"\n  weather = {wlab}:")
        print(f"    {'slice':12s} {'R2':>6s} {'[95% CI]':>16s}  "
              f"{'b_wx[z]':>12s}  {'b_gain[z]':>12s}")
        for s in C.SLICE_ORDER:
            key = f"v_{s}"
            if key not in sub.index:
                continue
            r = sub.loc[key]
            print(f"    {s:12s} {r.r2:6.3f} [{r.r2_lo:.3f}, {r.r2_hi:.3f}]  "
                  f"{r.beta_weather:+.3f}[{r.z_weather:+4.1f}]"
                  f"  {r.beta_course:+.3f}[{r.z_course:+4.1f}]  (n={int(r.n)})")
        present = [f"v_{s}" for s in C.SLICE_ORDER if f"v_{s}" in sub.index]
        if present:
            print(f"    -> highest R2: {sub.loc[present, 'r2'].idxmax()}")

    # ---- (B) simple metrics, ranked by R2 ----
    print("\n" + "=" * 86)
    print("B. Simple metrics ~ weather + total_gain  (ranked by R2, top 12)")
    print("=" * 86)
    met = out[out.kind == "metric"]
    for wlab in [w[1] for w in C.WEATHER_SETS]:
        sub = met[met.weather_set == wlab].sort_values("r2", ascending=False).head(12)
        print(f"\n  weather = {wlab}:")
        with pd.option_context("display.float_format", lambda v: f"{v:+.3f}"):
            print(sub[["target", "n", "r2", "beta_weather", "beta_course"]]
                  .to_string(index=False))

    # ---- (C) one combined ranking: slice v_j vs the best naive metric ----
    print("\n" + "=" * 86)
    print("C. Best-explained target overall (temp set), slices flagged with *")
    print("=" * 86)
    sub = out[out.weather_set == "temp"].sort_values("r2", ascending=False).head(15)
    for _, r in sub.iterrows():
        flag = " *" if r.kind == "slice" else "  "
        print(f"  {flag} {r.target:16s} R2={r.r2:.3f}  (n={int(r.n)})")

    print(f"\n[write] {pq}")


if __name__ == "__main__":
    main()
