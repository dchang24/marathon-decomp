"""q08 - which model's v_j best tracks the external conditions? (ALL_B)

Presentation step 2: holding the slice fixed at ALL_B, compare the four nested
production models by how strongly their race factor v_j correlates with the
objective race conditions (the q02 winner air temp + total elevation gain):

    baseline = u_i + v_j                         (rank-1 only)
    aging    = baseline + spline aging + gamma   (no per-athlete drift)
    drift    = baseline + d_i                     (no aging block)
    full     = baseline + aging + d_i             (production)

Each model's v_j is reloaded (no refit) from its registered ALL_B nu=8 fit and
re-expressed under the same G1 beta=0 APC gauge as the rest of the covariate
analysis. For each model we report:

  (A) univariate Pearson + Spearman (Fisher 95% CIs) of v_j vs air temp and vs
      total gain;
  (B) the joint 2-predictor standardized OLS  v_j ~ z(temp) + z(total_gain):
      R2 (with a 95% CI from the model's OWN bootstrap replicates -> athlete
      sampling noise) and the partial standardized betas (+ z).

If `full` shows the strongest covariate association, that supports it being the
better de-biased race factor (more of v_j is real physical difficulty, less is
field-composition / progression artifact). Console + results/analysis/covariate/
08_model_comparison/model_comparison__ALL_B.{parquet,csv}.

Run after q01 (needs covariate_merged for the covariates)::

    python scripts/05_analysis/covariate/q08_model_comparison.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

import covariate_common as C

SUBDIR = "08_model_comparison"
SLICE = "ALL_B"                 # the headline slice; models compared on it
WEATHER = "temp_field"          # q02 winner (air temp, C)
COURSE = C.COURSE_PREDICTOR     # total_gain_m


def fisher_ci(r: float, n: int, kind: str = "pearson",
              alpha: float = 0.05) -> tuple[float, float]:
    """95% CI for a correlation via the Fisher z-transform (Spearman SE x1.06)."""
    if not np.isfinite(r) or n < 4 or abs(r) >= 1.0:
        return np.nan, np.nan
    z = np.arctanh(r)
    se = (1.06 if kind == "spearman" else 1.0) / np.sqrt(n - 3)
    zc = stats.norm.ppf(1 - alpha / 2)
    return float(np.tanh(z - zc * se)), float(np.tanh(z + zc * se))


def boot_r2_ci(model: str, temp_by_rid: pd.Series, gain_by_rid: pd.Series,
               ) -> tuple[float, float, int]:
    """95% CI for the joint R2, from the model's own (re-gauged) replicates.

    Each replicate's v_j is regressed on the (fixed) standardized temp+gain over
    the joint complete-case races; we take the across-replicate R2 percentiles.
    Returns (lo, hi, R); (nan, nan, 0) if the model has no bootstrap.
    """
    res = C.load_vj_boot_wide(SLICE, model=model)
    if res is None:
        return np.nan, np.nan, 0
    _fdir, race_ids, _t, Bg = res                 # Bg: (R, J) re-gauged replicates
    temp = temp_by_rid.reindex(race_ids).to_numpy(float)
    gain = gain_by_rid.reindex(race_ids).to_numpy(float)
    finite_rep = np.isfinite(Bg).all(axis=0)      # races present in every replicate
    cc = np.isfinite(temp) & np.isfinite(gain) & finite_rep
    Xz = np.column_stack([C.standardize(temp[cc]), C.standardize(gain[cc])])
    X = np.column_stack([np.ones(int(cc.sum())), Xz])
    XtXinv = np.linalg.inv(X.T @ X)
    r2s = []
    for r in range(Bg.shape[0]):
        yz = C.standardize(Bg[r, cc])
        beta = XtXinv @ (X.T @ yz)
        resid = yz - X @ beta
        r2s.append(1 - (resid @ resid) / ((yz - yz.mean()) ** 2).sum())
    r2s = np.asarray(r2s)
    return float(np.percentile(r2s, 2.5)), float(np.percentile(r2s, 97.5)), Bg.shape[0]


def corr_pair(v: np.ndarray, x: np.ndarray, label: str) -> dict:
    m = np.isfinite(v) & np.isfinite(x)
    vv, xx = v[m], x[m]
    n = int(m.sum())
    pr, _ = stats.pearsonr(vv, xx)
    sr, _ = stats.spearmanr(vv, xx)
    pr_lo, pr_hi = fisher_ci(pr, n, "pearson")
    sr_lo, sr_hi = fisher_ci(sr, n, "spearman")
    return {f"pearson_{label}": pr, f"pearson_{label}_lo": pr_lo,
            f"pearson_{label}_hi": pr_hi, f"spearman_{label}": sr,
            f"spearman_{label}_lo": sr_lo, f"spearman_{label}_hi": sr_hi}


def main() -> None:
    merged = pd.read_parquet(C.MERGED_PATH)
    cov = merged.set_index("race_id")[[WEATHER, COURSE]]
    temp_by_rid, gain_by_rid = cov[WEATHER], cov[COURSE]

    rows = []
    for model in C.MODELS:
        d = C.load_slice_vj(SLICE, model=model).rename(columns={f"v_{SLICE}": "v"})
        d = d[["race_id", "v"]].merge(
            cov, left_on="race_id", right_index=True, how="left")
        v = d["v"].to_numpy(float)
        temp = d[WEATHER].to_numpy(float)
        gain = d[COURSE].to_numpy(float)

        # joint standardized OLS on the complete cases
        cc = np.isfinite(v) & np.isfinite(temp) & np.isfinite(gain)
        n = int(cc.sum())
        yz = C.standardize(v[cc])
        Xz = np.column_stack([C.standardize(temp[cc]), C.standardize(gain[cc])])
        X = np.column_stack([np.ones(n), Xz])
        beta, se, r2 = C.ols(X, yz)
        z = beta / se
        r2_lo, r2_hi, R = boot_r2_ci(model, temp_by_rid, gain_by_rid)

        row = dict(model=model, n=n, r2=r2, r2_lo=r2_lo, r2_hi=r2_hi, R_boot=R,
                   b_temp=beta[1], z_temp=z[1], b_gain=beta[2], z_gain=z[2])
        row.update(corr_pair(v, temp, "temp"))
        row.update(corr_pair(v, gain, "gain"))
        rows.append(row)

    out = pd.DataFrame(rows)
    pq = C.out_path(SUBDIR, "model_comparison", SLICE, "parquet")
    out.to_parquet(pq, index=False)
    out.to_csv(C.out_path(SUBDIR, "model_comparison", SLICE, "csv"), index=False)

    # -------------------- console report -------------------- #
    print("=" * 92)
    print(f"q08 model comparison - v_j vs temp + total_gain on {SLICE} "
          f"(full model, beta=0 gauge)")
    print("   does the production 'full' model's v_j track external conditions best?")
    print("=" * 92)
    print(f"\n  (A) univariate correlation of v_j with each covariate "
          f"(Spearman [95% CI]):")
    print(f"    {'model':9s} {'n':>4s}  {'rho(temp)':>20s}  {'rho(gain)':>20s}")
    for _, r in out.iterrows():
        ct = f"{r.spearman_temp:+.3f}[{r.spearman_temp_lo:+.2f},{r.spearman_temp_hi:+.2f}]"
        cg = f"{r.spearman_gain:+.3f}[{r.spearman_gain_lo:+.2f},{r.spearman_gain_hi:+.2f}]"
        print(f"    {r.model:9s} {int(r.n):4d}  {ct:>20s}  {cg:>20s}")

    print(f"\n  (B) joint OLS  v_j ~ z(temp) + z(gain)   "
          f"(R2 [95% CI from model bootstrap]; betas with z):")
    print(f"    {'model':9s} {'R2':>6s} {'[95% CI]':>16s}  "
          f"{'b_temp[z]':>13s}  {'b_gain[z]':>13s}")
    for _, r in out.iterrows():
        ci = (f"[{r.r2_lo:.3f},{r.r2_hi:.3f}]" if np.isfinite(r.r2_lo)
              else "[  -  ,  -  ]")
        print(f"    {r.model:9s} {r.r2:6.3f} {ci:>16s}  "
              f"{r.b_temp:+.3f}[{r.z_temp:+4.1f}]  {r.b_gain:+.3f}[{r.z_gain:+4.1f}]")

    best_r2 = out.loc[out.r2.idxmax(), "model"]
    best_rho = out.loc[out.spearman_temp.abs().idxmax(), "model"]
    print(f"\n  -> highest joint R2   : {best_r2}")
    print(f"  -> strongest |rho(temp)|: {best_rho}")
    print(f"\n[write] {pq}")


if __name__ == "__main__":
    main()
