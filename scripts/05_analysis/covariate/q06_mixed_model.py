"""q06 - random-course mixed model for v_j ~ temp + total_gain (ALL_B).

Method 2 of the "over-counting" fixes, and the principled one. A two-level model
with a **random intercept per physical course**:

    v_ij = b0 + b_temp * temp_ij + b_gain * gain_j + u_course + eps_ij
    u_course ~ N(0, tau^2),   eps_ij ~ N(0, sigma^2)

The random intercept absorbs the correlation among editions of the same course,
so partial pooling automatically discounts the pseudo-replicates and gives
`b_gain` (a course-level / level-2 covariate) a standard error with the correct
effective df. `b_temp` varies *within* course and is barely affected.

NOTE `total_gain_m` is constant within a course, so a course *fixed* effect would
be collinear with it (the effect would vanish). The course effect must be
**random** for `b_gain` to remain estimable -- this is exactly the right tool.

All variables standardized (betas comparable to q03-q05). Reports the fixed
effects with mixed-model SEs next to naive OLS, the variance components, and the
intraclass correlation. Console + results/analysis/covariate/06_mixed_model/
mixed_model__ALL_B.csv.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

import covariate_common as C

SUBDIR = "06_mixed_model"
SLICE = "ALL_B"             # this analysis is pinned to the ALL_B full fit


def main() -> None:
    df = pd.read_parquet(C.MERGED_PATH)
    y, X, cl, d = C.allb_design(df)
    n = len(y)
    G = pd.unique(cl).size

    dd = pd.DataFrame({"y": y, "temp": X[:, 1], "gain": X[:, 2], "course": cl})

    # naive OLS for side-by-side
    beta_ols = np.linalg.solve(X.T @ X, X.T @ y)
    resid = y - X @ beta_ols
    se_ols = np.sqrt(np.diag((resid @ resid) / (n - 3) * np.linalg.inv(X.T @ X)))

    # random-intercept mixed model (REML). `gain` is constant within course
    # (a pure level-2 covariate), which makes the REML surface flat near tau=0;
    # lbfgs gets stuck at the singular boundary, so use bfgs (cg/powell fall
    # back) and assert we left the boundary.
    md = smf.mixedlm("y ~ temp + gain", dd, groups=dd["course"])
    res = md.fit(reml=True, method=["bfgs", "cg", "powell"])

    tau2 = float(res.cov_re.iloc[0, 0])      # course-intercept variance
    if tau2 < 1e-6:
        raise SystemExit("mixed model degenerated to tau^2=0 (optimizer boundary); "
                         "inspect convergence before trusting SEs.")
    sig2 = float(res.scale)                  # residual variance
    icc = tau2 / (tau2 + sig2)

    print("=" * 84)
    print("q06 random-course mixed model - v_ALL_B ~ temp + gain + (1|course)")
    print(f"   n={n} editions, G={G} physical courses (random intercept), REML")
    print("=" * 84)
    print(f"\n  variance components: tau^2(course)={tau2:.4f}  "
          f"sigma^2(resid)={sig2:.4f}  ICC={icc:.3f}")
    print(f"  (ICC = {100 * icc:.0f}% of residual v_j variance is between-course)")

    print(f"\n  {'term':6s} {'OLS b[se]':>16s}   {'MixedLM b[se]':>18s}   "
          f"{'z':>5s}  {'p':>9s}   {'SE ratio':>8s}")
    rows = []
    for term, j in [("temp", 1), ("gain", 2)]:
        b_m = float(res.fe_params[term])
        se_m = float(res.bse_fe[term])
        z = b_m / se_m
        p = 2 * stats.norm.sf(abs(z))
        ratio = se_m / se_ols[j]
        print(f"  {term:6s} {beta_ols[j]:+.3f}[{se_ols[j]:.3f}]   "
              f"{b_m:+.3f}[{se_m:.3f}]      {z:+5.1f}  {p:9.1e}   x{ratio:5.2f}")
        rows.append(dict(term=term, b_ols=beta_ols[j], se_ols=se_ols[j],
                         b_mixed=b_m, se_mixed=se_m, z=z, p=p, se_ratio=ratio,
                         ci_lo=b_m - 1.96 * se_m, ci_hi=b_m + 1.96 * se_m))

    print(f"\n  -> gain SE inflates x{rows[1]['se_ratio']:.1f} vs iid OLS; "
          f"z {beta_ols[2] / se_ols[2]:+.1f} -> {rows[1]['z']:+.1f}")

    out = pd.DataFrame(rows)
    out.attrs = {}
    meta = pd.DataFrame([dict(term="_var", b_ols=np.nan, se_ols=np.nan,
                              b_mixed=tau2, se_mixed=sig2, z=icc, p=np.nan,
                              se_ratio=np.nan, ci_lo=np.nan, ci_hi=np.nan)])
    dst = C.out_path(SUBDIR, "mixed_model", SLICE, "csv")
    pd.concat([out, meta], ignore_index=True).to_csv(dst, index=False)
    print(f"\n[write] {dst}")


if __name__ == "__main__":
    main()
