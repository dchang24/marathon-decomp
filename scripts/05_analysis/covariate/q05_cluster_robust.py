"""q05 - cluster-robust inference for v_j ~ temp + total_gain (ALL_B).

Method 1 of the "over-counting" fixes. `total_gain_m` varies only at the
**physical-course** level (series_key x course_id): ALL_B has 343 editions but
only ~62 distinct courses, so editions on the same course are pseudo-replicates
of one gain value. The naive iid OLS SE treats all 343 as independent and
massively over-states the precision of `b_gain` (it is fine for `b_temp`, which
varies within course).

We refit the **same standardized OLS** (point estimate unchanged) and re-price
the SEs four ways, clustering on the physical course:

    iid        naive OLS                         (q03 baseline, over-confident)
    CR1        cluster sandwich, G/(G-1) corr    (Liang-Zeger)
    CR2        Bell-McCaffrey bias-reduced       (recommended for ~62 clusters)
    pairs-boot resample whole courses w/ repl    (the unit is the course)

t-CIs for CR1/CR2 use dof = G-1 (conservative). The WLS point fit is also shown
with iid vs CR1 SEs, to check whether its stronger `b_gain` survives honest
inference. Console + results/analysis/covariate/05_cluster_robust/
cluster_robust__ALL_B.csv.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

import covariate_common as C

SUBDIR = "05_cluster_robust"
SLICE = "ALL_B"             # this analysis is pinned to the ALL_B full fit
PREDS = ["intercept", "temp", "gain"]


def ols_beta(X, y):
    return np.linalg.solve(X.T @ X, X.T @ y)


def cov_iid(X, resid):
    n, k = X.shape
    sigma2 = (resid @ resid) / (n - k)
    return sigma2 * np.linalg.inv(X.T @ X)


def cov_cluster(X, resid, groups, kind="CR1"):
    """Cluster-robust sandwich. kind in {CR0, CR1, CR2}."""
    n, k = X.shape
    bread = np.linalg.inv(X.T @ X)
    meat = np.zeros((k, k))
    for g in groups:
        Xg, ug = X[g], resid[g]
        if kind == "CR2":
            Hgg = Xg @ bread @ Xg.T
            A = np.eye(len(g)) - Hgg
            w, Q = np.linalg.eigh(A)
            w = np.clip(w, 1e-10, None)
            ug = (Q @ np.diag(w ** -0.5) @ Q.T) @ ug
        s = Xg.T @ ug
        meat += np.outer(s, s)
    V = bread @ meat @ bread
    if kind == "CR1":
        G = len(groups)
        V *= G / (G - 1) * (n - 1) / (n - k)
    return V


def pairs_cluster_boot(X, y, groups, *, B=4000, seed=0):
    """Resample whole courses with replacement -> (B, k) beta draws."""
    rng = np.random.default_rng(seed)
    G = len(groups)
    out = []
    for _ in range(B):
        pick = rng.integers(0, G, G)
        idx = np.concatenate([groups[p] for p in pick])
        Xb, yb = X[idx], y[idx]
        try:
            out.append(np.linalg.solve(Xb.T @ Xb, Xb.T @ yb))
        except np.linalg.LinAlgError:
            continue
    return np.array(out)


def report(name, beta, V, dof, fh):
    se = np.sqrt(np.diag(V))
    tc = stats.t.ppf(0.975, dof)
    for j in (1, 2):
        z = beta[j] / se[j]
        p = 2 * stats.t.sf(abs(z), dof)
        lo, hi = beta[j] - tc * se[j], beta[j] + tc * se[j]
        print(f"    {PREDS[j]:5s}  b={beta[j]:+.3f}  se={se[j]:.3f}  "
              f"t={z:+5.1f}  p={p:7.1e}  95%CI[{lo:+.3f},{hi:+.3f}]")
        fh.append(dict(fit=name, term=PREDS[j], beta=beta[j], se=se[j],
                       t=z, p=p, ci_lo=lo, ci_hi=hi, dof=dof))


def main() -> None:
    df = pd.read_parquet(C.MERGED_PATH)
    y, X, cl, d = C.allb_design(df)
    n = len(y)
    uniq = pd.unique(cl)
    groups = [np.where(cl == u)[0] for u in uniq]
    G = len(groups)
    sizes = np.array([len(g) for g in groups])
    print("=" * 80)
    print("q05 cluster-robust inference - v_ALL_B ~ temp + total_gain (standardized)")
    print(f"   n={n} editions in G={G} physical courses "
          f"(median {int(np.median(sizes))}/course, max {sizes.max()})")
    print("=" * 80)

    beta = ols_beta(X, y)
    resid = y - X @ beta
    rows: list[dict] = []

    print("\n  OLS point estimate, SE priced four ways:")
    print("\n  [iid]  naive OLS (treats 343 editions as independent):")
    report("iid", beta, cov_iid(X, resid), n - 3, rows)
    print("\n  [CR1]  cluster-robust (Liang-Zeger), dof=G-1:")
    report("CR1", beta, cov_cluster(X, resid, groups, "CR1"), G - 1, rows)
    print("\n  [CR2]  Bell-McCaffrey bias-reduced, dof=G-1:")
    report("CR2", beta, cov_cluster(X, resid, groups, "CR2"), G - 1, rows)

    print("\n  [pairs-boot]  resample whole courses (B=4000):")
    bb = pairs_cluster_boot(X, y, groups)
    for j in (1, 2):
        lo, hi = np.percentile(bb[:, j], [2.5, 97.5])
        se = bb[:, j].std(ddof=1)
        print(f"    {PREDS[j]:5s}  b={beta[j]:+.3f}  se={se:.3f}  "
              f"95%CI[{lo:+.3f},{hi:+.3f}]  (sign-stable: "
              f"{100 * np.mean(np.sign(bb[:, j]) == np.sign(beta[j])):.0f}%)")
        rows.append(dict(fit="pairs_boot", term=PREDS[j], beta=beta[j], se=se,
                         t=np.nan, p=np.nan, ci_lo=lo, ci_hi=hi, dof=G - 1))

    # ---- WLS point fit, iid vs cluster-robust ---- #
    sd = C.load_vj_boot_sd("ALL_B")
    d2 = d.merge(sd, on="race_id", how="left")
    var = d2["vsd_ALL_B"].to_numpy(float) ** 2
    pw = 1.0 / np.maximum(var, 1e-12)
    pw = pw / pw.mean()
    Xw = X * np.sqrt(pw)[:, None]
    yw = y * np.sqrt(pw)
    beta_w = ols_beta(Xw, yw)
    resid_w = yw - Xw @ beta_w
    print("\n  WLS (1/var_j) point estimate, iid vs cluster-robust:")
    print("    [WLS-iid]:")
    report("WLS_iid", beta_w, cov_iid(Xw, resid_w), n - 3, rows)
    print("    [WLS-CR2]:")
    report("WLS_CR2", beta_w, cov_cluster(Xw, resid_w, groups, "CR2"), G - 1, rows)

    out = C.out_path(SUBDIR, "cluster_robust", SLICE, "csv")
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"\n[write] {out}")


if __name__ == "__main__":
    main()
