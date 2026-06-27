# Covariate validation: v_j vs weather + course (ALL_B)

*Does the race factor v_j track objective race-day conditions? Post-hoc on the fitted v_j (full production model, beta=0 gauge); no refit.*

*Single-slice analyses (q02/q05/q06/q08) are pinned to ALL_B; the q03/q04 blocks highlight the chosen slice.*

_Generated 2026-06-21._


## VARIABLE SELECTION  v_j vs each covariate  (q02, ALL_B)

| Quantity | Value | Source |
|---|---|---|
| weather air temp (C) | Pearson +0.653  Spearman +0.634 [+0.56,+0.70]  (n=343)  <= best | `02_variable_selection/variable_selection__ALL_B.csv, row temp_field` |
| weather WBGT field (C) | Pearson +0.628  Spearman +0.624 [+0.55,+0.69]  (n=343) | `02_variable_selection/variable_selection__ALL_B.csv, row wbgt_field` |
| weather WBGT max (C) | Pearson +0.529  Spearman +0.528 [+0.44,+0.60]  (n=343) | `02_variable_selection/variable_selection__ALL_B.csv, row wbgt_max` |
| weather precip (mm) | Pearson -0.017  Spearman -0.127 [-0.24,-0.02]  (n=343) | `02_variable_selection/variable_selection__ALL_B.csv, row precip_total_mm` |
| weather wind (m/s) | Pearson -0.087  Spearman -0.095 [-0.21,+0.02]  (n=343) | `02_variable_selection/variable_selection__ALL_B.csv, row wind_field` |
| weather WBGT rising slope (C/h) | Pearson +0.038  Spearman +0.084 [-0.03,+0.19]  (n=343) | `02_variable_selection/variable_selection__ALL_B.csv, row wbgt_rising_slope_c_per_h` |
| course  total elevation gain (m) | Pearson +0.272  Spearman +0.332 [+0.23,+0.43]  (n=343)  <= best | `02_variable_selection/variable_selection__ALL_B.csv, row total_gain_m` |
| course  \|net elevation gain\| (m) | Pearson +0.203  Spearman +0.299 [+0.19,+0.40]  (n=343) | `02_variable_selection/variable_selection__ALL_B.csv, row abs_net_gain_m` |
| course  net elevation gain (m) | Pearson -0.158  Spearman -0.186 [-0.29,-0.08]  (n=343) | `02_variable_selection/variable_selection__ALL_B.csv, row net_gain_m` |

## JOINT REGRESSION  v_j ~ temp+gain  (q03)

| Quantity | Value | Source |
|---|---|---|
| ALL_B R2 [95% CI] | 0.458 [0.383,0.535]  b_temp +0.627(z+15.5)  b_gain +0.177(z+4.4)  <= highlight | `03_regression/regression__6slices.csv, weather_set=temp, target=v_ALL_B` |
| ALL_M R2 [95% CI] | 0.448 [0.374,0.521]  b_temp +0.617(z+15.1)  b_gain +0.182(z+4.5) | `03_regression/regression__6slices.csv, weather_set=temp, target=v_ALL_M` |
| ALL_W R2 [95% CI] | 0.396 [0.310,0.479]  b_temp +0.593(z+13.8)  b_gain +0.140(z+3.3) | `03_regression/regression__6slices.csv, weather_set=temp, target=v_ALL_W` |
| Po10_B R2 [95% CI] | 0.493 [0.391,0.584]  b_temp +0.653(z+14.8)  b_gain +0.208(z+4.7) | `03_regression/regression__6slices.csv, weather_set=temp, target=v_Po10_B` |
| Po10_M R2 [95% CI] | 0.511 [0.412,0.607]  b_temp +0.660(z+14.8)  b_gain +0.228(z+5.1) | `03_regression/regression__6slices.csv, weather_set=temp, target=v_Po10_M` |
| Po10_W R2 [95% CI] | 0.483 [0.357,0.594]  b_temp +0.623(z+12.2)  b_gain +0.269(z+5.3) | `03_regression/regression__6slices.csv, weather_set=temp, target=v_Po10_W` |

## JOINT REGRESSION  v_j ~ WBGT+gain  (q03)

| Quantity | Value | Source |
|---|---|---|
| ALL_B R2 [95% CI] | 0.331 [0.265,0.408]  b_WBGT +0.509(z+11.4)  b_gain +0.225(z+5.1)  <= highlight | `03_regression/regression__6slices.csv, weather_set=WBGT, target=v_ALL_B` |
| ALL_M R2 [95% CI] | 0.331 [0.267,0.402]  b_WBGT +0.507(z+11.4)  b_gain +0.228(z+5.1) | `03_regression/regression__6slices.csv, weather_set=WBGT, target=v_ALL_M` |
| ALL_W R2 [95% CI] | 0.273 [0.197,0.359]  b_WBGT +0.472(z+10.1)  b_gain +0.185(z+4.0) | `03_regression/regression__6slices.csv, weather_set=WBGT, target=v_ALL_W` |
| Po10_B R2 [95% CI] | 0.281 [0.197,0.371]  b_WBGT +0.460(z+8.8)  b_gain +0.254(z+4.9) | `03_regression/regression__6slices.csv, weather_set=WBGT, target=v_Po10_B` |
| Po10_M R2 [95% CI] | 0.297 [0.208,0.397]  b_WBGT +0.468(z+8.8)  b_gain +0.277(z+5.2) | `03_regression/regression__6slices.csv, weather_set=WBGT, target=v_Po10_M` |
| Po10_W R2 [95% CI] | 0.276 [0.174,0.391]  b_WBGT +0.423(z+7.0)  b_gain +0.306(z+5.1) | `03_regression/regression__6slices.csv, weather_set=WBGT, target=v_Po10_W` |

## v_j vs SIMPLE METRICS  (joint ~ temp+gain; directly comparable)

| Quantity | Value | Source |
|---|---|---|
| v_ALL_B (full model)  [REF] | R2 0.458 [0.383,0.535]  b_temp +0.627(z+15.5)  b_gain +0.177(z+4.4) | `03_regression/regression__6slices.csv, target=v_ALL_B` |
| p1_time_sec | R2 0.141 [0.094,0.210]  b_temp +0.146(z+2.9)  b_gain +0.325(z+6.4) | `03_regression/regression__6slices.csv, kind=metric, target=p1_time_sec` |
| p10_time_sec | R2 0.093 [0.052,0.164]  b_temp +0.186(z+3.6)  b_gain +0.215(z+4.1) | `03_regression/regression__6slices.csv, kind=metric, target=p10_time_sec` |
| sub4_pct | R2 0.090 [0.040,0.163]  b_temp -0.252(z-4.8)  b_gain -0.130(z-2.5) | `03_regression/regression__6slices.csv, kind=metric, target=sub4_pct` |
| m_top25_sec | R2 0.085 [0.047,0.141]  b_temp +0.122(z+2.3)  b_gain +0.247(z+4.7) | `03_regression/regression__6slices.csv, kind=metric, target=m_top25_sec` |
| m_sub4_pct | R2 0.077 [0.022,0.170]  b_temp -0.252(z-4.8)  b_gain -0.084(z-1.6) | `03_regression/regression__6slices.csv, kind=metric, target=m_sub4_pct` |
| bq_pct | R2 0.076 [0.041,0.133]  b_temp -0.186(z-3.5)  b_gain -0.178(z-3.4) | `03_regression/regression__6slices.csv, kind=metric, target=bq_pct` |
| p25_time_sec | R2 0.063 [0.030,0.124]  b_temp +0.168(z+3.2)  b_gain +0.162(z+3.1) | `03_regression/regression__6slices.csv, kind=metric, target=p25_time_sec` |
| m_bq_pct | R2 0.063 [0.026,0.126]  b_temp -0.168(z-3.2)  b_gain -0.162(z-3.0) | `03_regression/regression__6slices.csv, kind=metric, target=m_bq_pct` |
| mean_time_sec | R2 0.057 [0.022,0.112]  b_temp +0.183(z+3.4)  b_gain +0.129(z+2.4) | `03_regression/regression__6slices.csv, kind=metric, target=mean_time_sec` |
| w_top1_sec | R2 0.056 [0.027,0.102]  b_temp -0.049(z-0.9)  b_gain +0.239(z+4.5) | `03_regression/regression__6slices.csv, kind=metric, target=w_top1_sec` |
| -> R2 ratio  v_ALL_B / best metric | 0.458 / 0.141 = 3.2x  (p1_time_sec) | `derived` |
| -> R2 ratio  v_ALL_B / mean_time_sec | 0.458 / 0.057 = 8.0x | `derived` |
| -> R2 ratio  v_ALL_B / median_time_sec | 0.458 / 0.047 = 9.7x | `derived` |

## MODEL COMPARISON  v_j ~ temp+gain  (q08, ALL_B)

| Quantity | Value | Source |
|---|---|---|
| baseline  R2 [95% CI] | 0.376 [0.346,0.387]  rho(temp) 0.551 | `08_model_comparison/model_comparison__ALL_B.csv, model=baseline` |
| aging     R2 [95% CI] | 0.409 [0.375,0.419]  rho(temp) 0.588 | `08_model_comparison/model_comparison__ALL_B.csv, model=aging` |
| drift     R2 [95% CI] | 0.443 [0.418,0.457]  rho(temp) 0.609 | `08_model_comparison/model_comparison__ALL_B.csv, model=drift` |
| full      R2 [95% CI] | 0.458 [0.432,0.471]  rho(temp) 0.634 | `08_model_comparison/model_comparison__ALL_B.csv, model=full` |
| -> aging gain (baseline->aging) | +0.033 | `derived` |
| -> d_i gain  (baseline->drift) | +0.067 | `derived` |
| -> d_i gain  (aging->full) | +0.049 | `derived` |
| -> best model | full (R2 0.458) | `derived` |

## ROBUSTNESS  v_j ~ temp+gain, 2x2 fits  (q04, ALL_B)

| Quantity | Value | Source |
|---|---|---|
| OLS | b_temp +0.627(z+15.5)  b_gain +0.177(z+4.4)  R2 0.458 | `04_regression_robustness/regression_robustness__6slices.csv, target=v_ALL_B, fit=OLS` |
| WLS | b_temp +0.672(z+16.9)  b_gain +0.419(z+10.2)  R2 0.552 | `04_regression_robustness/regression_robustness__6slices.csv, target=v_ALL_B, fit=WLS` |
| Robust | b_temp +0.592(z+16.9)  b_gain +0.158(z+4.5)  R2 0.497 | `04_regression_robustness/regression_robustness__6slices.csv, target=v_ALL_B, fit=Robust` |
| Robust+WLS | b_temp +0.624(z+18.1)  b_gain +0.470(z+12.7)  R2 0.598 | `04_regression_robustness/regression_robustness__6slices.csv, target=v_ALL_B, fit=Robust+WLS` |

## OVER-COUNTING #1: cluster-robust SE  (q05, ALL_B)

| Quantity | Value | Source |
|---|---|---|
| physical courses G (clusters) | 62  (vs n=343 editions) | `05_cluster_robust/cluster_robust__ALL_B.csv` |
| temp z: iid -> CR2 | 15.5 -> 11.9  (robust) | `05_cluster_robust/cluster_robust__ALL_B.csv` |
| gain z: iid -> CR1 -> CR2 | 4.4 -> 1.7 -> 1.6 | `05_cluster_robust/cluster_robust__ALL_B.csv` |
| gain z (WLS): iid -> CR2 | 10.2 -> 1.8 | `05_cluster_robust/cluster_robust__ALL_B.csv` |

## OVER-COUNTING #2: random-course mixed model  (q06, ALL_B)

| Quantity | Value | Source |
|---|---|---|
| temp: OLS z -> mixed z | +15.5 -> +18.2  (SE x0.91) | `06_mixed_model/mixed_model__ALL_B.csv, term=temp` |
| gain: OLS z -> mixed z | +4.4 -> +2.2  (SE x1.71) | `06_mixed_model/mixed_model__ALL_B.csv, term=gain` |
| ICC (between-course share of resid var) | 0.365 | `06_mixed_model/mixed_model__ALL_B.csv, _var row (z=ICC)` |
