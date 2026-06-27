"""q02 - which weather / course variable best tracks race difficulty?

Variable selection against the single ALL_B **full**-model v_j (the headline
fit, per request). For every candidate covariate we report Pearson + Spearman
(with Fisher 95% CIs) vs gauged v_j, ranked by |Spearman|.

  weather candidates: air temp, WBGT field, WBGT max, wind, precip, WBGT slope
  course  candidates: total gain, net gain, |net gain|

Higher v_j = harder race, so a positive correlation = "more of this covariate =
slower race". The legacy pass found air temp ~ WBGT both dominant (air temp
marginally ahead); this confirms / refutes that on production data.

Console + results/analysis/covariate/02_variable_selection/
variable_selection__ALL_B.{parquet,csv}.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

import covariate_common as C

SUBDIR = "02_variable_selection"
TARGET = f"v_{C.VAR_SELECT_SLICE}"


def fisher_ci(r: float, n: int, kind: str = "pearson",
              alpha: float = 0.05) -> tuple[float, float]:
    """95% CI for a correlation via the Fisher z-transform.

    Spearman uses the Fieller-Hartley-Pearson SE inflation 1.06 / sqrt(n-3).
    """
    if not np.isfinite(r) or n < 4 or abs(r) >= 1.0:
        return np.nan, np.nan
    z = np.arctanh(r)
    se = (1.06 if kind == "spearman" else 1.0) / np.sqrt(n - 3)
    zc = stats.norm.ppf(1 - alpha / 2)
    return float(np.tanh(z - zc * se)), float(np.tanh(z + zc * se))


def corr_row(x, y, cov_col, cov_lab, family) -> dict:
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = len(x)
    if n < 5 or np.std(x) == 0 or np.std(y) == 0:
        return dict(covariate=cov_col, label=cov_lab, family=family, n=n,
                    pearson=np.nan, p_pear=np.nan, pearson_lo=np.nan, pearson_hi=np.nan,
                    spearman=np.nan, p_spear=np.nan, spearman_lo=np.nan, spearman_hi=np.nan)
    pr, pp = stats.pearsonr(x, y)
    sr, sp = stats.spearmanr(x, y)
    pr_lo, pr_hi = fisher_ci(pr, n, "pearson")
    sr_lo, sr_hi = fisher_ci(sr, n, "spearman")
    return dict(covariate=cov_col, label=cov_lab, family=family, n=n,
                pearson=pr, p_pear=pp, pearson_lo=pr_lo, pearson_hi=pr_hi,
                spearman=sr, p_spear=sp, spearman_lo=sr_lo, spearman_hi=sr_hi)


def main() -> None:
    df = pd.read_parquet(C.MERGED_PATH)
    y = df[TARGET].to_numpy(float)

    rows = []
    for col, lab in C.WEATHER_CANDIDATES:
        rows.append(corr_row(df[col].to_numpy(float), y, col, lab, "weather"))
    for col, lab in C.COURSE_CANDIDATES:
        rows.append(corr_row(df[col].to_numpy(float), y, col, lab, "course"))
    res = pd.DataFrame(rows)
    res["abs_spearman"] = res.spearman.abs()
    pq = C.out_path(SUBDIR, "variable_selection", C.VAR_SELECT_SLICE, "parquet")
    res.to_parquet(pq, index=False)
    res.to_csv(C.out_path(SUBDIR, "variable_selection", C.VAR_SELECT_SLICE, "csv"),
               index=False)

    print("=" * 80)
    print(f"q02 variable selection - covariate vs {TARGET} (full model, beta=0 gauge)")
    print("   higher v_j = harder race; ranked by |Spearman| within family")
    print("=" * 80)
    for family in ("weather", "course"):
        sub = res[res.family == family].sort_values("abs_spearman", ascending=False)
        print(f"\n  {family}:")
        print(f"    {'covariate':24s} {'n':>4s}  {'Pearson':>8s}   "
              f"{'Spearman [95% CI]':>26s}")
        for _, r in sub.iterrows():
            ci = (f"[{r.spearman_lo:+.2f}, {r.spearman_hi:+.2f}]"
                  if np.isfinite(r.spearman_lo) else "[   nan,    nan]")
            print(f"    {r.label:24s} {int(r.n):4d}  {r.pearson:+8.3f}   "
                  f"{r.spearman:+.3f} {ci}")
        best = sub.iloc[0]
        print(f"    -> best {family}: {best.label}  "
              f"(|Spearman|={best.abs_spearman:.3f})")

    # head-to-head on the two legacy front-runners
    wx = res.set_index("covariate")
    if {"temp_field", "wbgt_max"}.issubset(wx.index):
        t, w = wx.loc["temp_field"], wx.loc["wbgt_max"]
        print("\n" + "-" * 80)
        print("  air temp vs WBGT max (legacy: both best, air temp marginally ahead)")
        print(f"    air temp  Pearson {t.pearson:+.3f}  Spearman {t.spearman:+.3f}")
        print(f"    WBGT max  Pearson {w.pearson:+.3f}  Spearman {w.spearman:+.3f}")
        lead = "air temp" if t.abs_spearman >= abs(w.spearman) else "WBGT max"
        print(f"    -> stronger |Spearman|: {lead}")

    print(f"\n[write] {pq}")


if __name__ == "__main__":
    main()
