# Baseline noise model: Student-t selection (ALL_B_14-25_mrc2)

*Rank-1 baseline log t = u_i + v_j; ALS/Anderson share the fixed point.*

*Degrees of freedom nu fixed by held-out predictive density (the in-sample fit always prefers heavier tails); a single shared nu=8 is used.*

_Generated 2026-06-21._


## WHY STUDENT-T

| Quantity | Value | Source |
|---|---|---|
| sigma^2  (L2, nu=inf) | 3.1693e-03 | `ALL_B_14-25_mrc2/param_sensitivity.csv (solver=anderson), nu=inf, col sigma2` |
| sigma^2  (Student-t, nu=8) | 1.8833e-03 | `ALL_B_14-25_mrc2/param_sensitivity.csv (solver=anderson), nu=8, col sigma2` |
| -> L2 variance inflation | 1.683x  (68% larger) | `derived` |
| QQ R^2  baseline L2 vs Normal | 0.9235 | `qq_plot/residual_diagnostics.csv, model='baseline (u+v), L2', qq_r2_normal` |
| QQ R^2  baseline t8 vs t(8) | 0.9792 | `qq_plot/residual_diagnostics.csv, model='baseline (u+v), t8', qq_r2_ref` |
| skewness (baseline t8 resid) | 1.356 | `qq_plot/residual_diagnostics.csv, model='baseline (u+v), t8', skewness` |
| excess kurtosis (baseline t8) | 11.940 | `qq_plot/residual_diagnostics.csv, model='baseline (u+v), t8', excess_kurtosis` |
| nu implied by kurtosis | 4.503 | `qq_plot/residual_diagnostics.csv, nu_implied_by_kurtosis  (=6/exk+4)` |

## IN-SAMPLE FIT & COMPLEXITY (AIC / BIC / eff. d.o.f.)

| Quantity | Value | Source |
|---|---|---|
| log-lik   (L2, nu=inf) | 1,814,393.8 | `ALL_B_14-25_mrc2/param_sensitivity.csv (solver=anderson), nu=inf, full_loglik` |
| log-lik   (Student-t, nu=8) | 1,914,942.6 | `ALL_B_14-25_mrc2/param_sensitivity.csv (solver=anderson), nu=8, full_loglik` |
| effective d.o.f.  (L2 / t8) | 402,943 / 402,943  (=AIC/2+loglik; equal, nu not counted) | `derived` |
| AIC       (L2 / t8) | -2,822,902 / -3,023,999 | `ALL_B_14-25_mrc2/param_sensitivity.csv (solver=anderson), full_aic` |
| -> dAIC (t8 - L2) | -201,098  (t8 lower = better) | `derived` |
| BIC       (L2 / t8) | 2,026,145 / 1,825,047 | `ALL_B_14-25_mrc2/param_sensitivity.csv (solver=anderson), full_bic` |
| -> dBIC (t8 - L2) | -201,098  (t8 lower; dof unchanged) | `derived` |

## CHOOSING NU (out-of-sample CV)

| Quantity | Value | Source |
|---|---|---|
| held-out logdens/cell  (L2, nu=inf) | 0.8351 | `ALL_B_14-25_mrc2/nu_selection.csv (source=grid), nu=inf, cv_per_cell` |
| held-out logdens/cell  (nu=8) | 0.9077 | `ALL_B_14-25_mrc2/nu_selection.csv (source=grid), nu=8, cv_per_cell` |
| -> CV gain (nats/cell) | +0.0726 | `derived` |
| => term earns its keep? | in-sample BIC -201,098 AND held-out CV +0.0726  (both improve -> not in-sample overfit) | `derived` |
| grid argmax nu | 8 | `ALL_B_14-25_mrc2/nu_selection.csv (source=grid), argmax cv_per_cell` |
| grid-best nu | 8 | `ALL_B_14-25_mrc2/selected_nu.csv, grid_best_nu` |
| Brent-refined nu* | 7.976 | `ALL_B_14-25_mrc2/selected_nu.csv, brent_nu` |
| selected (shared) nu | 8 | `ALL_B_14-25_mrc2/selected_nu.csv, selected_nu` |
| 1-SE acceptable nu interval | [8, 10] | `nu_decision.csv, this slice, nu_1se_lo/hi` |
| recommended shared nu (across slices) | 8 | `nu_decision.csv, slice=__RECOMMENDED__, nu_argmax` |
| v ranking stability vs L2 @nu=8 | Spearman 0.9909, Pearson 0.9926 | `ALL_B_14-25_mrc2/param_sensitivity.csv (solver=anderson), nu=8, spearman_vs_L2 / pearson_vs_L2` |
| max \|dv\| vs L2 (log-time) | 0.0174 | `ALL_B_14-25_mrc2/param_sensitivity.csv (solver=anderson), nu=8, max_abs_dv` |
