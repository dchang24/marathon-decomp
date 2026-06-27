# Aging vs Drift -- 2x2 term grid

2x2 aging-vs-drift term grid   nu=8   mrc=[2, 5]   slices requested: Po10_W, Po10_M, Po10_B, ALL_W, ALL_M, ALL_B, WA_M, WA_W

## Po10_W_14-25_mrc2   nu=8

I=17,679  J=203  N=61,907   cells: baseline / aging / drift / full

### 1. Fit 2x2 (descriptive: loglik / RSS / AIC / BIC)

| cell | label | loglik | rss | aic | bic | n_params |
| --- | --- | --- | --- | --- | --- | --- |
| baseline | no-d / no-aging | 93716.5 | 224.313 | -151671 | 9855.08 | 17882 |
| aging | no-d / +aging | 95768.3 | 213.475 | -155763 | 5817.68 | 17888 |
| drift | +d / no-aging | 104871 | 158.567 | -161190 | 58107.8 | 27010 |
| full | +d / +aging | 106024 | 156.287 | -163813 | 54053.9 | 27016 |


loglik on the 2x2 grid:
| loglik | no-aging | +aging |
| --- | --- | --- |
| no-d | 93716.4687 | 95768.2729 |
| +d | 104871.2643 | 106024.3980 |


marginal loglik gains (does a term survive adding the other?):
| effect | d_loglik |
| --- | --- |
| aging_gain \| no-d | 2051.8041 |
| aging_gain \| +d | 1153.1338 |
| d_gain \| no-aging | 11154.7955 |
| d_gain \| +aging | 10256.1252 |
| interaction (non-additivity) | -898.6704 |


### 2. Contribution overlap in the full fit (aging vs d_i)

| quantity | value |
| --- | --- |
| corr_aging_drift | -0.0767618 |
| var_aging | 0.000621991 |
| var_drift | 0.00071902 |
| cov | -5.13352e-05 |
| sd_ratio_d_over_aging | 1.07517 |


(corr near 0 => terms carry distinct signal => both needed; gauge-sensitive in the linear direction.)

### 3. Drift block: does turning aging ON shrink d_i?

| quantity | value |
| --- | --- |
| omega_d2_drift | 0.000327307 |
| omega_d2_full | 0.00027776 |
| n_eligible | 9127 |
| d_sd_drift | 0.015082 |
| d_sd_full | 0.0136845 |
| frac_improver_drift | 0.542347 |
| frac_improver_full | 0.557796 |
| corr_d_drift_full | 0.935368 |
| max_abs_dd | 0.0340855 |


### 4. Aging block: does turning d ON change the aging curve / gamma?

| quantity | value |
| --- | --- |
| curve_corr | 0.995944 |
| curve_max_abs_delta | 0.0304439 |
| peak_age_aging | 2.46332 |
| peak_age_full | 2.34602 |
| gamma_max_abs_delta | 0.000220927 |
| gamma_corr | 0.998581 |


### 5. Factor stability across the grid (v_j and u_i)


v_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9730 | 0.9251 | 0.9199 |
| aging | 0.9730 | 1.0000 | 0.9022 | 0.9202 |
| drift | 0.9251 | 0.9022 | 1.0000 | 0.9294 |
| full | 0.9199 | 0.9202 | 0.9294 | 1.0000 |

v_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9574 | 0.9060 | 0.8992 |
| aging | 0.9574 | 1.0000 | 0.8622 | 0.8745 |
| drift | 0.9060 | 0.8622 | 1.0000 | 0.8899 |
| full | 0.8992 | 0.8745 | 0.8899 | 1.0000 |

v_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.01784 | 0.02928 | 0.03532 |
| aging | 0.01784 | 0 | 0.02772 | 0.02437 |
| drift | 0.02928 | 0.02772 | 0 | 0.03199 |
| full | 0.03532 | 0.02437 | 0.03199 | 0 |


u_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9956 | 0.9962 | 0.9923 |
| aging | 0.9956 | 1.0000 | 0.9922 | 0.9961 |
| drift | 0.9962 | 0.9922 | 1.0000 | 0.9914 |
| full | 0.9923 | 0.9961 | 0.9914 | 1.0000 |

u_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9953 | 0.9967 | 0.9921 |
| aging | 0.9953 | 1.0000 | 0.9926 | 0.9963 |
| drift | 0.9967 | 0.9926 | 1.0000 | 0.9920 |
| full | 0.9921 | 0.9963 | 0.9920 | 1.0000 |

u_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.3354 | 0.4539 | 0.5041 |
| aging | 0.3354 | 0 | 0.4631 | 0.5154 |
| drift | 0.4539 | 0.4631 | 0 | 0.509 |
| full | 0.5041 | 0.5154 | 0.509 | 0 |


## Po10_W_14-25_mrc5   nu=8

I=3,457  J=166  N=25,316   cells: baseline / aging / drift / full

### 1. Fit 2x2 (descriptive: loglik / RSS / AIC / BIC)

| cell | label | loglik | rss | aic | bic | n_params |
| --- | --- | --- | --- | --- | --- | --- |
| baseline | no-d / no-aging | 35590.5 | 110.348 | -63937 | -34456.8 | 3623 |
| aging | no-d / +aging | 36808 | 102.504 | -66360 | -36831.1 | 3629 |
| drift | +d / no-aging | 41457.4 | 71.5882 | -69975.4 | -17317.5 | 7081 |
| full | +d / +aging | 41593.6 | 71.6898 | -70579.8 | -19272.4 | 7087 |


loglik on the 2x2 grid:
| loglik | no-aging | +aging |
| --- | --- | --- |
| no-d | 35590.5013 | 36808.0235 |
| +d | 41457.3556 | 41593.6386 |


marginal loglik gains (does a term survive adding the other?):
| effect | d_loglik |
| --- | --- |
| aging_gain \| no-d | 1217.5222 |
| aging_gain \| +d | 136.2831 |
| d_gain \| no-aging | 5866.8543 |
| d_gain \| +aging | 4785.6151 |
| interaction (non-additivity) | -1081.2391 |


### 2. Contribution overlap in the full fit (aging vs d_i)

| quantity | value |
| --- | --- |
| corr_aging_drift | -0.0542875 |
| var_aging | 0.00123804 |
| var_drift | 0.000877472 |
| cov | -5.6585e-05 |
| sd_ratio_d_over_aging | 0.841877 |


(corr near 0 => terms carry distinct signal => both needed; gauge-sensitive in the linear direction.)

### 3. Drift block: does turning aging ON shrink d_i?

| quantity | value |
| --- | --- |
| omega_d2_drift | 0.000278929 |
| omega_d2_full | 0.000189672 |
| n_eligible | 3457 |
| d_sd_drift | 0.0151602 |
| d_sd_full | 0.0121182 |
| frac_improver_drift | 0.491177 |
| frac_improver_full | 0.50188 |
| corr_d_drift_full | 0.873569 |
| max_abs_dd | 0.0275947 |


### 4. Aging block: does turning d ON change the aging curve / gamma?

| quantity | value |
| --- | --- |
| curve_corr | 0.999041 |
| curve_max_abs_delta | 0.0135637 |
| peak_age_aging | 2.58062 |
| peak_age_full | 2.69792 |
| gamma_max_abs_delta | 3.18647e-05 |
| gamma_corr | 0.999854 |


### 5. Factor stability across the grid (v_j and u_i)


v_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9259 | 0.8194 | 0.9028 |
| aging | 0.9259 | 1.0000 | 0.8222 | 0.9592 |
| drift | 0.8194 | 0.8222 | 1.0000 | 0.8675 |
| full | 0.9028 | 0.9592 | 0.8675 | 1.0000 |

v_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.8948 | 0.7889 | 0.8666 |
| aging | 0.8948 | 1.0000 | 0.7882 | 0.9386 |
| drift | 0.7889 | 0.7882 | 1.0000 | 0.8311 |
| full | 0.8666 | 0.9386 | 0.8311 | 1.0000 |

v_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.02196 | 0.05009 | 0.02603 |
| aging | 0.02196 | 0 | 0.04535 | 0.02228 |
| drift | 0.05009 | 0.04535 | 0 | 0.04165 |
| full | 0.02603 | 0.02228 | 0.04165 | 0 |


u_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9867 | 0.9973 | 0.9837 |
| aging | 0.9867 | 1.0000 | 0.9864 | 0.9980 |
| drift | 0.9973 | 0.9864 | 1.0000 | 0.9869 |
| full | 0.9837 | 0.9980 | 0.9869 | 1.0000 |

u_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9817 | 0.9969 | 0.9778 |
| aging | 0.9817 | 1.0000 | 0.9821 | 0.9979 |
| drift | 0.9969 | 0.9821 | 1.0000 | 0.9814 |
| full | 0.9778 | 0.9979 | 0.9814 | 1.0000 |

u_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.133 | 0.1636 | 0.1464 |
| aging | 0.133 | 0 | 0.1721 | 0.1549 |
| drift | 0.1636 | 0.1721 | 0 | 0.1299 |
| full | 0.1464 | 0.1549 | 0.1299 | 0 |


### 6. Aging-curve de-biasing across mrc -- Po10_W (mrc 2 vs 5)

population aging curve, centered to mean 0 on shared A_n in [0, 11.67] yr (100 pts), one APC beta=0 gauge.

peak age (A_n at fastest, yr) per curve:
| mrc | variant | peak_age |
| --- | --- | --- |
| 2 | noD | 2.476 |
| 2 | withD | 2.358 |
| 5 | noD | 2.594 |
| 5 | withD | 2.712 |


curve-distance contrasts (RMS on the shared grid):
| contrast | rms |
| --- | --- |
| gap_noD (everyone vs dedicated, no d) | 0.0051946 |
| gap_withD (everyone vs dedicated, +d) | 0.0023468 |
| d_on_mrc2 (\|\|aging-full\|\| @ everyone) | 0.011023 |
| d_on_mrc5 (\|\|aging-full\|\| @ dedicated) | 0.0053158 |


(expect gap_withD << gap_noD and d_on_mrc2 >> d_on_mrc5 if no-d/everyone is the contaminated curve.)

## Po10_M_14-25_mrc2   nu=8

I=32,023  J=252  N=121,520   cells: baseline / aging / drift / full

### 1. Fit 2x2 (descriptive: loglik / RSS / AIC / BIC)

| cell | label | loglik | rss | aic | bic | n_params |
| --- | --- | --- | --- | --- | --- | --- |
| baseline | no-d / no-aging | 164759 | 618.166 | -264971 | 48339.7 | 32275 |
| aging | no-d / +aging | 169520 | 583.11 | -274479 | 38889.7 | 32281 |
| drift | +d / no-aging | 186076 | 447.877 | -282816 | 150809 | 50434 |
| full | +d / +aging | 188468 | 441.155 | -288648 | 139898 | 50440 |


loglik on the 2x2 grid:
| loglik | no-aging | +aging |
| --- | --- | --- |
| no-d | 164759.4567 | 169519.5697 |
| +d | 186075.6573 | 188468.3335 |


marginal loglik gains (does a term survive adding the other?):
| effect | d_loglik |
| --- | --- |
| aging_gain \| no-d | 4760.1130 |
| aging_gain \| +d | 2392.6762 |
| d_gain \| no-aging | 21316.2006 |
| d_gain \| +aging | 18948.7638 |
| interaction (non-additivity) | -2367.4368 |


### 2. Contribution overlap in the full fit (aging vs d_i)

| quantity | value |
| --- | --- |
| corr_aging_drift | -0.076981 |
| var_aging | 0.000994248 |
| var_drift | 0.000892283 |
| cov | -7.2508e-05 |
| sd_ratio_d_over_aging | 0.947336 |


(corr near 0 => terms carry distinct signal => both needed; gauge-sensitive in the linear direction.)

### 3. Drift block: does turning aging ON shrink d_i?

| quantity | value |
| --- | --- |
| omega_d2_drift | 0.000363015 |
| omega_d2_full | 0.000286736 |
| n_eligible | 18158 |
| d_sd_drift | 0.0157008 |
| d_sd_full | 0.0136374 |
| frac_improver_drift | 0.539817 |
| frac_improver_full | 0.542075 |
| corr_d_drift_full | 0.912328 |
| max_abs_dd | 0.0514256 |


### 4. Aging block: does turning d ON change the aging curve / gamma?

| quantity | value |
| --- | --- |
| curve_corr | 0.998898 |
| curve_max_abs_delta | 0.0309836 |
| peak_age_aging | 2.64361 |
| peak_age_full | 2.64361 |
| gamma_max_abs_delta | 2.61152e-05 |
| gamma_corr | 0.999865 |


### 5. Factor stability across the grid (v_j and u_i)


v_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9681 | 0.9391 | 0.9439 |
| aging | 0.9681 | 1.0000 | 0.9001 | 0.9379 |
| drift | 0.9391 | 0.9001 | 1.0000 | 0.9220 |
| full | 0.9439 | 0.9379 | 0.9220 | 1.0000 |

v_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9575 | 0.9229 | 0.9378 |
| aging | 0.9575 | 1.0000 | 0.8743 | 0.9162 |
| drift | 0.9229 | 0.8743 | 1.0000 | 0.8944 |
| full | 0.9378 | 0.9162 | 0.8944 | 1.0000 |

v_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.01888 | 0.03362 | 0.03836 |
| aging | 0.01888 | 0 | 0.03969 | 0.03494 |
| drift | 0.03362 | 0.03969 | 0 | 0.03938 |
| full | 0.03836 | 0.03494 | 0.03938 | 0 |


u_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9933 | 0.9962 | 0.9899 |
| aging | 0.9933 | 1.0000 | 0.9901 | 0.9959 |
| drift | 0.9962 | 0.9901 | 1.0000 | 0.9889 |
| full | 0.9899 | 0.9959 | 0.9889 | 1.0000 |

u_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9929 | 0.9967 | 0.9899 |
| aging | 0.9929 | 1.0000 | 0.9906 | 0.9964 |
| drift | 0.9967 | 0.9906 | 1.0000 | 0.9900 |
| full | 0.9899 | 0.9964 | 0.9900 | 1.0000 |

u_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.3965 | 0.4418 | 0.5677 |
| aging | 0.3965 | 0 | 0.4526 | 0.513 |
| drift | 0.4418 | 0.4526 | 0 | 0.5802 |
| full | 0.5677 | 0.513 | 0.5802 | 0 |


## Po10_M_14-25_mrc5   nu=8

I=7,592  J=216  N=57,728   cells: baseline / aging / drift / full

### 1. Fit 2x2 (descriptive: loglik / RSS / AIC / BIC)

| cell | label | loglik | rss | aic | bic | n_params |
| --- | --- | --- | --- | --- | --- | --- |
| baseline | no-d / no-aging | 73637.5 | 338.577 | -131661 | -61683 | 7808 |
| aging | no-d / +aging | 76716.7 | 310.401 | -137807 | -67775.7 | 7814 |
| drift | +d / no-aging | 86478 | 226.072 | -144915 | -19242.1 | 15401 |
| full | +d / +aging | 86880.1 | 225.544 | -146640 | -25095 | 15407 |


loglik on the 2x2 grid:
| loglik | no-aging | +aging |
| --- | --- | --- |
| no-d | 73637.4907 | 76716.7322 |
| +d | 86478.0289 | 86880.1083 |


marginal loglik gains (does a term survive adding the other?):
| effect | d_loglik |
| --- | --- |
| aging_gain \| no-d | 3079.2415 |
| aging_gain \| +d | 402.0795 |
| d_gain \| no-aging | 12840.5382 |
| d_gain \| +aging | 10163.3762 |
| interaction (non-additivity) | -2677.1620 |


### 2. Contribution overlap in the full fit (aging vs d_i)

| quantity | value |
| --- | --- |
| corr_aging_drift | -0.047713 |
| var_aging | 0.00174589 |
| var_drift | 0.00104456 |
| cov | -6.44345e-05 |
| sd_ratio_d_over_aging | 0.773497 |


(corr near 0 => terms carry distinct signal => both needed; gauge-sensitive in the linear direction.)

### 3. Drift block: does turning aging ON shrink d_i?

| quantity | value |
| --- | --- |
| omega_d2_drift | 0.000326835 |
| omega_d2_full | 0.000205599 |
| n_eligible | 7592 |
| d_sd_drift | 0.0163563 |
| d_sd_full | 0.0124762 |
| frac_improver_drift | 0.494468 |
| frac_improver_full | 0.50382 |
| corr_d_drift_full | 0.845147 |
| max_abs_dd | 0.0330777 |


### 4. Aging block: does turning d ON change the aging curve / gamma?

| quantity | value |
| --- | --- |
| curve_corr | 0.998249 |
| curve_max_abs_delta | 0.0142193 |
| peak_age_aging | 2.64361 |
| peak_age_full | 2.8786 |
| gamma_max_abs_delta | 8.78677e-05 |
| gamma_corr | 0.997961 |


### 5. Factor stability across the grid (v_j and u_i)


v_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.8970 | 0.8483 | 0.8994 |
| aging | 0.8970 | 1.0000 | 0.8157 | 0.9709 |
| drift | 0.8483 | 0.8157 | 1.0000 | 0.8477 |
| full | 0.8994 | 0.9709 | 0.8477 | 1.0000 |

v_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.8587 | 0.8276 | 0.8591 |
| aging | 0.8587 | 1.0000 | 0.8000 | 0.9687 |
| drift | 0.8276 | 0.8000 | 1.0000 | 0.8352 |
| full | 0.8591 | 0.9687 | 0.8352 | 1.0000 |

v_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.036 | 0.05816 | 0.0433 |
| aging | 0.036 | 0 | 0.0555 | 0.02816 |
| drift | 0.05816 | 0.0555 | 0 | 0.05013 |
| full | 0.0433 | 0.02816 | 0.05013 | 0 |


u_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9844 | 0.9971 | 0.9813 |
| aging | 0.9844 | 1.0000 | 0.9848 | 0.9981 |
| drift | 0.9971 | 0.9848 | 1.0000 | 0.9848 |
| full | 0.9813 | 0.9981 | 0.9848 | 1.0000 |

u_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9775 | 0.9963 | 0.9737 |
| aging | 0.9775 | 1.0000 | 0.9790 | 0.9978 |
| drift | 0.9963 | 0.9790 | 1.0000 | 0.9786 |
| full | 0.9737 | 0.9978 | 0.9786 | 1.0000 |

u_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.1587 | 0.1967 | 0.1896 |
| aging | 0.1587 | 0 | 0.2164 | 0.1908 |
| drift | 0.1967 | 0.2164 | 0 | 0.2371 |
| full | 0.1896 | 0.1908 | 0.2371 | 0 |


### 6. Aging-curve de-biasing across mrc -- Po10_M (mrc 2 vs 5)

population aging curve, centered to mean 0 on shared A_n in [0, 11.69] yr (100 pts), one APC beta=0 gauge.

peak age (A_n at fastest, yr) per curve:
| mrc | variant | peak_age |
| --- | --- | --- |
| 2 | noD | 2.598 |
| 2 | withD | 2.716 |
| 5 | noD | 2.598 |
| 5 | withD | 2.952 |


curve-distance contrasts (RMS on the shared grid):
| contrast | rms |
| --- | --- |
| gap_noD (everyone vs dedicated, no d) | 0.0055027 |
| gap_withD (everyone vs dedicated, +d) | 0.0025684 |
| d_on_mrc2 (\|\|aging-full\|\| @ everyone) | 0.01135 |
| d_on_mrc5 (\|\|aging-full\|\| @ dedicated) | 0.0058791 |


(expect gap_withD << gap_noD and d_on_mrc2 >> d_on_mrc5 if no-d/everyone is the contaminated curve.)

## Po10_B_14-25_mrc2   nu=8

I=49,871  J=265  N=184,448   cells: baseline / aging / drift / full

### 1. Fit 2x2 (descriptive: loglik / RSS / AIC / BIC)

| cell | label | loglik | rss | aic | bic | n_params |
| --- | --- | --- | --- | --- | --- | --- |
| baseline | no-d / no-aging | 259357 | 848.047 | -418444 | 89179.5 | 50136 |
| aging | no-d / +aging | 266174 | 802.391 | -432065 | 75618.3 | 50142 |
| drift | +d / no-aging | 292011 | 611.241 | -446020 | 252625 | 77555 |
| full | +d / +aging | 295580 | 601.687 | -454492 | 237397 | 77561 |


loglik on the 2x2 grid:
| loglik | no-aging | +aging |
| --- | --- | --- |
| no-d | 259356.7917 | 266173.7284 |
| +d | 292011.3906 | 295579.7624 |


marginal loglik gains (does a term survive adding the other?):
| effect | d_loglik |
| --- | --- |
| aging_gain \| no-d | 6816.9368 |
| aging_gain \| +d | 3568.3717 |
| d_gain \| no-aging | 32654.5990 |
| d_gain \| +aging | 29406.0339 |
| interaction (non-additivity) | -3248.5650 |


### 2. Contribution overlap in the full fit (aging vs d_i)

| quantity | value |
| --- | --- |
| corr_aging_drift | -0.0713115 |
| var_aging | 0.000833873 |
| var_drift | 0.000839669 |
| cov | -5.96714e-05 |
| sd_ratio_d_over_aging | 1.00347 |


(corr near 0 => terms carry distinct signal => both needed; gauge-sensitive in the linear direction.)

### 3. Drift block: does turning aging ON shrink d_i?

| quantity | value |
| --- | --- |
| omega_d2_drift | 0.000350019 |
| omega_d2_full | 0.000284757 |
| n_eligible | 27418 |
| d_sd_drift | 0.0154709 |
| d_sd_full | 0.0136815 |
| frac_improver_drift | 0.542016 |
| frac_improver_full | 0.548946 |
| corr_d_drift_full | 0.920807 |
| max_abs_dd | 0.0460848 |


### 4. Aging block: does turning d ON change the aging curve / gamma?

| quantity | value |
| --- | --- |
| curve_corr | 0.996875 |
| curve_max_abs_delta | 0.0313383 |
| peak_age_aging | 2.7611 |
| peak_age_full | 2.64361 |
| gamma_max_abs_delta | 4.11602e-05 |
| gamma_corr | 0.999998 |


### 5. Factor stability across the grid (v_j and u_i)


v_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9749 | 0.9405 | 0.9269 |
| aging | 0.9749 | 1.0000 | 0.9124 | 0.9290 |
| drift | 0.9405 | 0.9124 | 1.0000 | 0.9167 |
| full | 0.9269 | 0.9290 | 0.9167 | 1.0000 |

v_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9668 | 0.9234 | 0.9062 |
| aging | 0.9668 | 1.0000 | 0.8806 | 0.8918 |
| drift | 0.9234 | 0.8806 | 1.0000 | 0.8760 |
| full | 0.9062 | 0.8918 | 0.8760 | 1.0000 |

v_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.0166 | 0.03103 | 0.03957 |
| aging | 0.0166 | 0 | 0.03251 | 0.03218 |
| drift | 0.03103 | 0.03251 | 0 | 0.03972 |
| full | 0.03957 | 0.03218 | 0.03972 | 0 |


u_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9951 | 0.9963 | 0.9919 |
| aging | 0.9951 | 1.0000 | 0.9922 | 0.9963 |
| drift | 0.9963 | 0.9922 | 1.0000 | 0.9907 |
| full | 0.9919 | 0.9963 | 0.9907 | 1.0000 |

u_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9950 | 0.9967 | 0.9921 |
| aging | 0.9950 | 1.0000 | 0.9925 | 0.9967 |
| drift | 0.9967 | 0.9925 | 1.0000 | 0.9913 |
| full | 0.9921 | 0.9967 | 0.9913 | 1.0000 |

u_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.439 | 0.4881 | 0.5569 |
| aging | 0.439 | 0 | 0.504 | 0.5462 |
| drift | 0.4881 | 0.504 | 0 | 0.5642 |
| full | 0.5569 | 0.5462 | 0.5642 | 0 |


## Po10_B_14-25_mrc5   nu=8

I=11,248  J=236  N=84,876   cells: baseline / aging / drift / full

### 1. Fit 2x2 (descriptive: loglik / RSS / AIC / BIC)

| cell | label | loglik | rss | aic | bic | n_params |
| --- | --- | --- | --- | --- | --- | --- |
| baseline | no-d / no-aging | 111488 | 458.893 | -200009 | -92655.1 | 11484 |
| aging | no-d / +aging | 115871 | 422.219 | -208763 | -101353 | 11490 |
| drift | +d / no-aging | 130557 | 304.587 | -219720 | -26228.1 | 22733 |
| full | +d / +aging | 131093 | 304.027 | -222092 | -34672.9 | 22739 |


loglik on the 2x2 grid:
| loglik | no-aging | +aging |
| --- | --- | --- |
| no-d | 111487.5432 | 115870.5801 |
| +d | 130556.5457 | 131092.8537 |


marginal loglik gains (does a term survive adding the other?):
| effect | d_loglik |
| --- | --- |
| aging_gain \| no-d | 4383.0368 |
| aging_gain \| +d | 536.3080 |
| d_gain \| no-aging | 19069.0025 |
| d_gain \| +aging | 15222.2737 |
| interaction (non-additivity) | -3846.7289 |


### 2. Contribution overlap in the full fit (aging vs d_i)

| quantity | value |
| --- | --- |
| corr_aging_drift | -0.0481979 |
| var_aging | 0.00156419 |
| var_drift | 0.000991636 |
| cov | -6.00282e-05 |
| sd_ratio_d_over_aging | 0.796216 |


(corr near 0 => terms carry distinct signal => both needed; gauge-sensitive in the linear direction.)

### 3. Drift block: does turning aging ON shrink d_i?

| quantity | value |
| --- | --- |
| omega_d2_drift | 0.000312154 |
| omega_d2_full | 0.00020067 |
| n_eligible | 11248 |
| d_sd_drift | 0.0159912 |
| d_sd_full | 0.0123569 |
| frac_improver_drift | 0.494488 |
| frac_improver_full | 0.504623 |
| corr_d_drift_full | 0.853414 |
| max_abs_dd | 0.0319804 |


### 4. Aging block: does turning d ON change the aging curve / gamma?

| quantity | value |
| --- | --- |
| curve_corr | 0.998481 |
| curve_max_abs_delta | 0.0136836 |
| peak_age_aging | 2.64361 |
| peak_age_full | 2.8786 |
| gamma_max_abs_delta | 7.34907e-05 |
| gamma_corr | 0.998129 |


### 5. Factor stability across the grid (v_j and u_i)


v_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9123 | 0.8501 | 0.9019 |
| aging | 0.9123 | 1.0000 | 0.8284 | 0.9625 |
| drift | 0.8501 | 0.8284 | 1.0000 | 0.8457 |
| full | 0.9019 | 0.9625 | 0.8457 | 1.0000 |

v_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.8755 | 0.8372 | 0.8822 |
| aging | 0.8755 | 1.0000 | 0.8162 | 0.9630 |
| drift | 0.8372 | 0.8162 | 1.0000 | 0.8419 |
| full | 0.8822 | 0.9630 | 0.8419 | 1.0000 |

v_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.03366 | 0.05558 | 0.05447 |
| aging | 0.03366 | 0 | 0.05204 | 0.05129 |
| drift | 0.05558 | 0.05204 | 0 | 0.04858 |
| full | 0.05447 | 0.05129 | 0.04858 | 0 |


u_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9872 | 0.9975 | 0.9846 |
| aging | 0.9872 | 1.0000 | 0.9872 | 0.9983 |
| drift | 0.9975 | 0.9872 | 1.0000 | 0.9874 |
| full | 0.9846 | 0.9983 | 0.9874 | 1.0000 |

u_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9848 | 0.9973 | 0.9820 |
| aging | 0.9848 | 1.0000 | 0.9852 | 0.9982 |
| drift | 0.9973 | 0.9852 | 1.0000 | 0.9851 |
| full | 0.9820 | 0.9982 | 0.9851 | 1.0000 |

u_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.1514 | 0.2018 | 0.1975 |
| aging | 0.1514 | 0 | 0.2195 | 0.1962 |
| drift | 0.2018 | 0.2195 | 0 | 0.244 |
| full | 0.1975 | 0.1962 | 0.244 | 0 |


### 6. Aging-curve de-biasing across mrc -- Po10_B (mrc 2 vs 5)

population aging curve, centered to mean 0 on shared A_n in [0, 11.69] yr (100 pts), one APC beta=0 gauge.

peak age (A_n at fastest, yr) per curve:
| mrc | variant | peak_age |
| --- | --- | --- |
| 2 | noD | 2.716 |
| 2 | withD | 2.598 |
| 5 | noD | 2.598 |
| 5 | withD | 2.834 |


curve-distance contrasts (RMS on the shared grid):
| contrast | rms |
| --- | --- |
| gap_noD (everyone vs dedicated, no d) | 0.0062365 |
| gap_withD (everyone vs dedicated, +d) | 0.0022516 |
| d_on_mrc2 (\|\|aging-full\|\| @ everyone) | 0.011406 |
| d_on_mrc5 (\|\|aging-full\|\| @ dedicated) | 0.0056468 |


(expect gap_withD << gap_noD and d_on_mrc2 >> d_on_mrc5 if no-d/everyone is the contaminated curve.)

## ALL_W_14-25_mrc2   nu=8

I=124,650  J=337  N=381,031   cells: baseline / aging / drift / full

### 1. Fit 2x2 (descriptive: loglik / RSS / AIC / BIC)

| cell | label | loglik | rss | aic | bic | n_params |
| --- | --- | --- | --- | --- | --- | --- |
| baseline | no-d / no-aging | 601485 | 1176.99 | -952999 | 403179 | 124987 |
| aging | no-d / +aging | 608936 | 1232.37 | -967887 | 388355 | 124993 |
| drift | +d / no-aging | 662343 | 847.692 | -1.0036e+06 | 738381 | 176062 |
| full | +d / +aging | 665864 | 926.213 | -1.01208e+06 | 722085 | 176068 |


loglik on the 2x2 grid:
| loglik | no-aging | +aging |
| --- | --- | --- |
| no-d | 601485.3372 | 608935.6039 |
| +d | 662342.8534 | 665864.1125 |


marginal loglik gains (does a term survive adding the other?):
| effect | d_loglik |
| --- | --- |
| aging_gain \| no-d | 7450.2668 |
| aging_gain \| +d | 3521.2590 |
| d_gain \| no-aging | 60857.5162 |
| d_gain \| +aging | 56928.5085 |
| interaction (non-additivity) | -3929.0077 |


### 2. Contribution overlap in the full fit (aging vs d_i)

| quantity | value |
| --- | --- |
| corr_aging_drift | -0.0285392 |
| var_aging | 0.000502185 |
| var_drift | 0.000565814 |
| cov | -1.52129e-05 |
| sd_ratio_d_over_aging | 1.06146 |


(corr near 0 => terms carry distinct signal => both needed; gauge-sensitive in the linear direction.)

### 3. Drift block: does turning aging ON shrink d_i?

| quantity | value |
| --- | --- |
| omega_d2_drift | 0.000295355 |
| omega_d2_full | 0.000263246 |
| n_eligible | 51074 |
| d_sd_drift | 0.0143393 |
| d_sd_full | 0.0133785 |
| frac_improver_drift | 0.50231 |
| frac_improver_full | 0.529917 |
| corr_d_drift_full | 0.964524 |
| max_abs_dd | 0.0373921 |


### 4. Aging block: does turning d ON change the aging curve / gamma?

| quantity | value |
| --- | --- |
| curve_corr | 0.991265 |
| curve_max_abs_delta | 0.0233578 |
| peak_age_aging | 2.52197 |
| peak_age_full | 2.34602 |
| gamma_max_abs_delta | 4.87138e-05 |
| gamma_corr | 0.999997 |


### 5. Factor stability across the grid (v_j and u_i)


v_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9582 | 0.9577 | 0.9448 |
| aging | 0.9582 | 1.0000 | 0.9260 | 0.9654 |
| drift | 0.9577 | 0.9260 | 1.0000 | 0.9387 |
| full | 0.9448 | 0.9654 | 0.9387 | 1.0000 |

v_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9540 | 0.9441 | 0.9363 |
| aging | 0.9540 | 1.0000 | 0.9090 | 0.9524 |
| drift | 0.9441 | 0.9090 | 1.0000 | 0.9259 |
| full | 0.9363 | 0.9524 | 0.9259 | 1.0000 |

v_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.02518 | 0.02508 | 0.02704 |
| aging | 0.02518 | 0 | 0.03297 | 0.02012 |
| drift | 0.02508 | 0.03297 | 0 | 0.02743 |
| full | 0.02704 | 0.02012 | 0.02743 | 0 |


u_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9977 | 0.9977 | 0.9958 |
| aging | 0.9977 | 1.0000 | 0.9957 | 0.9977 |
| drift | 0.9977 | 0.9957 | 1.0000 | 0.9953 |
| full | 0.9958 | 0.9977 | 0.9953 | 1.0000 |

u_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9975 | 0.9976 | 0.9956 |
| aging | 0.9975 | 1.0000 | 0.9956 | 0.9976 |
| drift | 0.9976 | 0.9956 | 1.0000 | 0.9953 |
| full | 0.9956 | 0.9976 | 0.9953 | 1.0000 |

u_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.4555 | 0.3877 | 0.465 |
| aging | 0.4555 | 0 | 0.4643 | 0.4222 |
| drift | 0.3877 | 0.4643 | 0 | 0.4739 |
| full | 0.465 | 0.4222 | 0.4739 | 0 |


## ALL_W_14-25_mrc5   nu=8

I=20,267  J=295  N=131,398   cells: baseline / aging / drift / full

### 1. Fit 2x2 (descriptive: loglik / RSS / AIC / BIC)

| cell | label | loglik | rss | aic | bic | n_params |
| --- | --- | --- | --- | --- | --- | --- |
| baseline | no-d / no-aging | 181454 | 570.038 | -321785 | -120576 | 20562 |
| aging | no-d / +aging | 184659 | 587.85 | -328185 | -126916 | 20568 |
| drift | +d / no-aging | 209510 | 367.758 | -346242 | 9865.97 | 40830 |
| full | +d / +aging | 209663 | 409.331 | -347656 | 3024.71 | 40836 |


loglik on the 2x2 grid:
| loglik | no-aging | +aging |
| --- | --- | --- |
| no-d | 181453.7386 | 184659.2867 |
| +d | 209510.4997 | 209663.1844 |


marginal loglik gains (does a term survive adding the other?):
| effect | d_loglik |
| --- | --- |
| aging_gain \| no-d | 3205.5481 |
| aging_gain \| +d | 152.6846 |
| d_gain \| no-aging | 28056.7611 |
| d_gain \| +aging | 25003.8977 |
| interaction (non-additivity) | -3052.8635 |


### 2. Contribution overlap in the full fit (aging vs d_i)

| quantity | value |
| --- | --- |
| corr_aging_drift | -0.0189295 |
| var_aging | 0.00121834 |
| var_drift | 0.000940031 |
| cov | -2.0258e-05 |
| sd_ratio_d_over_aging | 0.878389 |


(corr near 0 => terms carry distinct signal => both needed; gauge-sensitive in the linear direction.)

### 3. Drift block: does turning aging ON shrink d_i?

| quantity | value |
| --- | --- |
| omega_d2_drift | 0.000247727 |
| omega_d2_full | 0.000203771 |
| n_eligible | 20267 |
| d_sd_drift | 0.0139099 |
| d_sd_full | 0.0123902 |
| frac_improver_drift | 0.496324 |
| frac_improver_full | 0.507327 |
| corr_d_drift_full | 0.939984 |
| max_abs_dd | 0.0254254 |


### 4. Aging block: does turning d ON change the aging curve / gamma?

| quantity | value |
| --- | --- |
| curve_corr | 0.998913 |
| curve_max_abs_delta | 0.0116693 |
| peak_age_aging | 2.40467 |
| peak_age_full | 2.46332 |
| gamma_max_abs_delta | 8.34074e-06 |
| gamma_corr | 0.999991 |


### 5. Factor stability across the grid (v_j and u_i)


v_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9502 | 0.9067 | 0.9426 |
| aging | 0.9502 | 1.0000 | 0.8549 | 0.9572 |
| drift | 0.9067 | 0.8549 | 1.0000 | 0.9158 |
| full | 0.9426 | 0.9572 | 0.9158 | 1.0000 |

v_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9317 | 0.8849 | 0.9261 |
| aging | 0.9317 | 1.0000 | 0.8108 | 0.9401 |
| drift | 0.8849 | 0.8108 | 1.0000 | 0.8869 |
| full | 0.9261 | 0.9401 | 0.8869 | 1.0000 |

v_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.02171 | 0.043 | 0.02377 |
| aging | 0.02171 | 0 | 0.03759 | 0.01904 |
| drift | 0.043 | 0.03759 | 0 | 0.02904 |
| full | 0.02377 | 0.01904 | 0.02904 | 0 |


u_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9934 | 0.9986 | 0.9922 |
| aging | 0.9934 | 1.0000 | 0.9925 | 0.9987 |
| drift | 0.9986 | 0.9925 | 1.0000 | 0.9936 |
| full | 0.9922 | 0.9987 | 0.9936 | 1.0000 |

u_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9946 | 0.9985 | 0.9933 |
| aging | 0.9946 | 1.0000 | 0.9937 | 0.9987 |
| drift | 0.9985 | 0.9937 | 1.0000 | 0.9948 |
| full | 0.9933 | 0.9987 | 0.9948 | 1.0000 |

u_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.891 | 0.1737 | 0.8986 |
| aging | 0.891 | 0 | 0.8899 | 0.1883 |
| drift | 0.1737 | 0.8899 | 0 | 0.8977 |
| full | 0.8986 | 0.1883 | 0.8977 | 0 |


### 6. Aging-curve de-biasing across mrc -- ALL_W (mrc 2 vs 5)

population aging curve, centered to mean 0 on shared A_n in [0, 11.67] yr (100 pts), one APC beta=0 gauge.

peak age (A_n at fastest, yr) per curve:
| mrc | variant | peak_age |
| --- | --- | --- |
| 2 | noD | 2.594 |
| 2 | withD | 2.358 |
| 5 | noD | 2.358 |
| 5 | withD | 2.476 |


curve-distance contrasts (RMS on the shared grid):
| contrast | rms |
| --- | --- |
| gap_noD (everyone vs dedicated, no d) | 0.0057394 |
| gap_withD (everyone vs dedicated, +d) | 0.001972 |
| d_on_mrc2 (\|\|aging-full\|\| @ everyone) | 0.0080843 |
| d_on_mrc5 (\|\|aging-full\|\| @ dedicated) | 0.004224 |


(expect gap_withD << gap_noD and d_on_mrc2 >> d_on_mrc5 if no-d/everyone is the contaminated curve.)

## ALL_M_14-25_mrc2   nu=8

I=277,903  J=344  N=863,072   cells: baseline / aging / drift / full

### 1. Fit 2x2 (descriptive: loglik / RSS / AIC / BIC)

| cell | label | loglik | rss | aic | bic | n_params |
| --- | --- | --- | --- | --- | --- | --- |
| baseline | no-d / no-aging | 1.31452e+06 | 3061.74 | -2.07255e+06 | 1.1741e+06 | 278247 |
| aging | no-d / +aging | 1.33686e+06 | 3331.66 | -2.11721e+06 | 1.1295e+06 | 278253 |
| drift | +d / no-aging | 1.45614e+06 | 2199.61 | -2.18942e+06 | 2.02787e+06 | 399270 |
| full | +d / +aging | 1.46468e+06 | 2525.35 | -2.2119e+06 | 1.9739e+06 | 399276 |


loglik on the 2x2 grid:
| loglik | no-aging | +aging |
| --- | --- | --- |
| no-d | 1314520.1910 | 1336858.6607 |
| +d | 1456141.2265 | 1464680.9462 |


marginal loglik gains (does a term survive adding the other?):
| effect | d_loglik |
| --- | --- |
| aging_gain \| no-d | 22338.4697 |
| aging_gain \| +d | 8539.7197 |
| d_gain \| no-aging | 141621.0355 |
| d_gain \| +aging | 127822.2855 |
| interaction (non-additivity) | -13798.7500 |


### 2. Contribution overlap in the full fit (aging vs d_i)

| quantity | value |
| --- | --- |
| corr_aging_drift | -0.0133425 |
| var_aging | 0.000771201 |
| var_drift | 0.000618061 |
| cov | -9.21167e-06 |
| sd_ratio_d_over_aging | 0.895224 |


(corr near 0 => terms carry distinct signal => both needed; gauge-sensitive in the linear direction.)

### 3. Drift block: does turning aging ON shrink d_i?

| quantity | value |
| --- | --- |
| omega_d2_drift | 0.000331578 |
| omega_d2_full | 0.000281492 |
| n_eligible | 121022 |
| d_sd_drift | 0.0150968 |
| d_sd_full | 0.0136646 |
| frac_improver_drift | 0.500306 |
| frac_improver_full | 0.526202 |
| corr_d_drift_full | 0.955087 |
| max_abs_dd | 0.0407819 |


### 4. Aging block: does turning d ON change the aging curve / gamma?

| quantity | value |
| --- | --- |
| curve_corr | 0.995983 |
| curve_max_abs_delta | 0.0202187 |
| peak_age_aging | 2.7611 |
| peak_age_full | 2.52611 |
| gamma_max_abs_delta | 6.78648e-05 |
| gamma_corr | 0.999995 |


### 5. Factor stability across the grid (v_j and u_i)


v_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9515 | 0.9541 | 0.9384 |
| aging | 0.9515 | 1.0000 | 0.9113 | 0.9698 |
| drift | 0.9541 | 0.9113 | 1.0000 | 0.9091 |
| full | 0.9384 | 0.9698 | 0.9091 | 1.0000 |

v_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9586 | 0.9385 | 0.9342 |
| aging | 0.9586 | 1.0000 | 0.9072 | 0.9529 |
| drift | 0.9385 | 0.9072 | 1.0000 | 0.9049 |
| full | 0.9342 | 0.9529 | 0.9049 | 1.0000 |

v_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.03447 | 0.02909 | 0.03081 |
| aging | 0.03447 | 0 | 0.04418 | 0.01835 |
| drift | 0.02909 | 0.04418 | 0 | 0.04124 |
| full | 0.03081 | 0.01835 | 0.04124 | 0 |


u_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9979 | 0.9976 | 0.9958 |
| aging | 0.9979 | 1.0000 | 0.9960 | 0.9977 |
| drift | 0.9976 | 0.9960 | 1.0000 | 0.9952 |
| full | 0.9958 | 0.9977 | 0.9952 | 1.0000 |

u_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9978 | 0.9975 | 0.9956 |
| aging | 0.9978 | 1.0000 | 0.9959 | 0.9976 |
| drift | 0.9975 | 0.9959 | 1.0000 | 0.9952 |
| full | 0.9956 | 0.9976 | 0.9952 | 1.0000 |

u_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.3225 | 0.5373 | 0.5661 |
| aging | 0.3225 | 0 | 0.5683 | 0.5564 |
| drift | 0.5373 | 0.5683 | 0 | 0.575 |
| full | 0.5661 | 0.5564 | 0.575 | 0 |


## ALL_M_14-25_mrc5   nu=8

I=45,876  J=328  N=299,892   cells: baseline / aging / drift / full

### 1. Fit 2x2 (descriptive: loglik / RSS / AIC / BIC)

| cell | label | loglik | rss | aic | bic | n_params |
| --- | --- | --- | --- | --- | --- | --- |
| baseline | no-d / no-aging | 388643 | 1570.59 | -684881 | -194612 | 46204 |
| aging | no-d / +aging | 399285 | 1499.45 | -706152 | -215820 | 46210 |
| drift | +d / no-aging | 452158 | 1023.91 | -740375 | 129435 | 92081 |
| full | +d / +aging | 453221 | 1048.25 | -746178 | 104115 | 92087 |


loglik on the 2x2 grid:
| loglik | no-aging | +aging |
| --- | --- | --- |
| no-d | 388643.2709 | 399284.9643 |
| +d | 452158.3731 | 453220.9622 |


marginal loglik gains (does a term survive adding the other?):
| effect | d_loglik |
| --- | --- |
| aging_gain \| no-d | 10641.6934 |
| aging_gain \| +d | 1062.5892 |
| d_gain \| no-aging | 63515.1021 |
| d_gain \| +aging | 53935.9979 |
| interaction (non-additivity) | -9579.1042 |


### 2. Contribution overlap in the full fit (aging vs d_i)

| quantity | value |
| --- | --- |
| corr_aging_drift | -0.0214932 |
| var_aging | 0.00117483 |
| var_drift | 0.00102577 |
| cov | -2.35947e-05 |
| sd_ratio_d_over_aging | 0.934408 |


(corr near 0 => terms carry distinct signal => both needed; gauge-sensitive in the linear direction.)

### 3. Drift block: does turning aging ON shrink d_i?

| quantity | value |
| --- | --- |
| omega_d2_drift | 0.000290092 |
| omega_d2_full | 0.00022014 |
| n_eligible | 45876 |
| d_sd_drift | 0.0150393 |
| d_sd_full | 0.0127587 |
| frac_improver_drift | 0.499041 |
| frac_improver_full | 0.505362 |
| corr_d_drift_full | 0.911465 |
| max_abs_dd | 0.0328085 |


### 4. Aging block: does turning d ON change the aging curve / gamma?

| quantity | value |
| --- | --- |
| curve_corr | 0.99957 |
| curve_max_abs_delta | 0.0100071 |
| peak_age_aging | 2.46737 |
| peak_age_full | 2.52611 |
| gamma_max_abs_delta | 5.84078e-05 |
| gamma_corr | 0.997819 |


### 5. Factor stability across the grid (v_j and u_i)


v_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9208 | 0.9052 | 0.9200 |
| aging | 0.9208 | 1.0000 | 0.8397 | 0.9630 |
| drift | 0.9052 | 0.8397 | 1.0000 | 0.8750 |
| full | 0.9200 | 0.9630 | 0.8750 | 1.0000 |

v_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9011 | 0.8913 | 0.9017 |
| aging | 0.9011 | 1.0000 | 0.8005 | 0.9473 |
| drift | 0.8913 | 0.8005 | 1.0000 | 0.8503 |
| full | 0.9017 | 0.9473 | 0.8503 | 1.0000 |

v_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.02408 | 0.03884 | 0.03685 |
| aging | 0.02408 | 0 | 0.0369 | 0.02016 |
| drift | 0.03884 | 0.0369 | 0 | 0.03958 |
| full | 0.03685 | 0.02016 | 0.03958 | 0 |


u_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9931 | 0.9985 | 0.9916 |
| aging | 0.9931 | 1.0000 | 0.9925 | 0.9988 |
| drift | 0.9985 | 0.9925 | 1.0000 | 0.9932 |
| full | 0.9916 | 0.9988 | 0.9932 | 1.0000 |

u_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9933 | 0.9984 | 0.9917 |
| aging | 0.9933 | 1.0000 | 0.9927 | 0.9988 |
| drift | 0.9984 | 0.9927 | 1.0000 | 0.9934 |
| full | 0.9917 | 0.9988 | 0.9934 | 1.0000 |

u_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 1.07 | 0.1994 | 1.096 |
| aging | 1.07 | 0 | 1.06 | 0.192 |
| drift | 0.1994 | 1.06 | 0 | 1.086 |
| full | 1.096 | 0.192 | 1.086 | 0 |


### 6. Aging-curve de-biasing across mrc -- ALL_M (mrc 2 vs 5)

population aging curve, centered to mean 0 on shared A_n in [0, 11.69] yr (100 pts), one APC beta=0 gauge.

peak age (A_n at fastest, yr) per curve:
| mrc | variant | peak_age |
| --- | --- | --- |
| 2 | noD | 2.716 |
| 2 | withD | 2.480 |
| 5 | noD | 2.480 |
| 5 | withD | 2.598 |


curve-distance contrasts (RMS on the shared grid):
| contrast | rms |
| --- | --- |
| gap_noD (everyone vs dedicated, no d) | 0.0085199 |
| gap_withD (everyone vs dedicated, +d) | 0.0051563 |
| d_on_mrc2 (\|\|aging-full\|\| @ everyone) | 0.0071381 |
| d_on_mrc5 (\|\|aging-full\|\| @ dedicated) | 0.0038223 |


(expect gap_withD << gap_noD and d_on_mrc2 >> d_on_mrc5 if no-d/everyone is the contaminated curve.)

## ALL_B_14-25_mrc2   nu=8

I=402,597  J=347  N=1,244,290   cells: baseline / aging / drift / full

### 1. Fit 2x2 (descriptive: loglik / RSS / AIC / BIC)

| cell | label | loglik | rss | aic | bic | n_params |
| --- | --- | --- | --- | --- | --- | --- |
| baseline | no-d / no-aging | 1.91494e+06 | 4245.69 | -3.024e+06 | 1.82505e+06 | 402944 |
| aging | no-d / +aging | 1.94451e+06 | 4543.91 | -3.08312e+06 | 1.766e+06 | 402950 |
| drift | +d / no-aging | 2.11747e+06 | 3051.84 | -3.19158e+06 | 3.08642e+06 | 575072 |
| full | +d / +aging | 2.12945e+06 | 3431.03 | -3.22227e+06 | 3.01515e+06 | 575078 |


loglik on the 2x2 grid:
| loglik | no-aging | +aging |
| --- | --- | --- |
| no-d | 1914943.0364 | 1944506.7218 |
| +d | 2117474.5547 | 2129447.9812 |


marginal loglik gains (does a term survive adding the other?):
| effect | d_loglik |
| --- | --- |
| aging_gain \| no-d | 29563.6854 |
| aging_gain \| +d | 11973.4265 |
| d_gain \| no-aging | 202531.5183 |
| d_gain \| +aging | 184941.2594 |
| interaction (non-additivity) | -17590.2589 |


### 2. Contribution overlap in the full fit (aging vs d_i)

| quantity | value |
| --- | --- |
| corr_aging_drift | -0.0172698 |
| var_aging | 0.000664859 |
| var_drift | 0.000603597 |
| cov | -1.09402e-05 |
| sd_ratio_d_over_aging | 0.952816 |


(corr near 0 => terms carry distinct signal => both needed; gauge-sensitive in the linear direction.)

### 3. Drift block: does turning aging ON shrink d_i?

| quantity | value |
| --- | --- |
| omega_d2_drift | 0.000320584 |
| omega_d2_full | 0.000276648 |
| n_eligible | 172127 |
| d_sd_drift | 0.0148713 |
| d_sd_full | 0.0135982 |
| frac_improver_drift | 0.501101 |
| frac_improver_full | 0.527971 |
| corr_d_drift_full | 0.958329 |
| max_abs_dd | 0.0576559 |


### 4. Aging block: does turning d ON change the aging curve / gamma?

| quantity | value |
| --- | --- |
| curve_corr | 0.994837 |
| curve_max_abs_delta | 0.0211924 |
| peak_age_aging | 2.70236 |
| peak_age_full | 2.46737 |
| gamma_max_abs_delta | 2.46798e-05 |
| gamma_corr | 0.999995 |


### 5. Factor stability across the grid (v_j and u_i)


v_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9542 | 0.9556 | 0.9421 |
| aging | 0.9542 | 1.0000 | 0.9163 | 0.9687 |
| drift | 0.9556 | 0.9163 | 1.0000 | 0.9166 |
| full | 0.9421 | 0.9687 | 0.9166 | 1.0000 |

v_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9584 | 0.9418 | 0.9348 |
| aging | 0.9584 | 1.0000 | 0.9098 | 0.9531 |
| drift | 0.9418 | 0.9098 | 1.0000 | 0.9077 |
| full | 0.9348 | 0.9531 | 0.9077 | 1.0000 |

v_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.03112 | 0.02753 | 0.02788 |
| aging | 0.03112 | 0 | 0.0398 | 0.0187 |
| drift | 0.02753 | 0.0398 | 0 | 0.03502 |
| full | 0.02788 | 0.0187 | 0.03502 | 0 |


u_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9983 | 0.9978 | 0.9965 |
| aging | 0.9983 | 1.0000 | 0.9965 | 0.9980 |
| drift | 0.9978 | 0.9965 | 1.0000 | 0.9960 |
| full | 0.9965 | 0.9980 | 0.9960 | 1.0000 |

u_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9982 | 0.9977 | 0.9963 |
| aging | 0.9982 | 1.0000 | 0.9963 | 0.9980 |
| drift | 0.9977 | 0.9963 | 1.0000 | 0.9959 |
| full | 0.9963 | 0.9980 | 0.9959 | 1.0000 |

u_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.3323 | 0.5295 | 1.069 |
| aging | 0.3323 | 0 | 0.5236 | 0.9272 |
| drift | 0.5295 | 0.5236 | 0 | 1.048 |
| full | 1.069 | 0.9272 | 1.048 | 0 |


## ALL_B_14-25_mrc5   nu=8

I=66,271  J=330  N=432,253   cells: baseline / aging / drift / full

### 1. Fit 2x2 (descriptive: loglik / RSS / AIC / BIC)

| cell | label | loglik | rss | aic | bic | n_params |
| --- | --- | --- | --- | --- | --- | --- |
| baseline | no-d / no-aging | 570759 | 2147.93 | -1.00832e+06 | -277266 | 66601 |
| aging | no-d / +aging | 584515 | 2107.24 | -1.03582e+06 | -304700 | 66607 |
| drift | +d / no-aging | 662396 | 1397.18 | -1.08825e+06 | 209947 | 132873 |
| full | +d / +aging | 663556 | 1478.84 | -1.09536e+06 | 176616 | 132879 |


loglik on the 2x2 grid:
| loglik | no-aging | +aging |
| --- | --- | --- |
| no-d | 570759.0888 | 584515.4985 |
| +d | 662395.5249 | 663556.1258 |


marginal loglik gains (does a term survive adding the other?):
| effect | d_loglik |
| --- | --- |
| aging_gain \| no-d | 13756.4097 |
| aging_gain \| +d | 1160.6008 |
| d_gain \| no-aging | 91636.4361 |
| d_gain \| +aging | 79040.6272 |
| interaction (non-additivity) | -12595.8089 |


### 2. Contribution overlap in the full fit (aging vs d_i)

| quantity | value |
| --- | --- |
| corr_aging_drift | -0.0202267 |
| var_aging | 0.00122987 |
| var_drift | 0.00100086 |
| cov | -2.2441e-05 |
| sd_ratio_d_over_aging | 0.902105 |


(corr near 0 => terms carry distinct signal => both needed; gauge-sensitive in the linear direction.)

### 3. Drift block: does turning aging ON shrink d_i?

| quantity | value |
| --- | --- |
| omega_d2_drift | 0.000276528 |
| omega_d2_full | 0.000215445 |
| n_eligible | 66271 |
| d_sd_drift | 0.0146833 |
| d_sd_full | 0.0126564 |
| frac_improver_drift | 0.498016 |
| frac_improver_full | 0.505485 |
| corr_d_drift_full | 0.920732 |
| max_abs_dd | 0.0336286 |


### 4. Aging block: does turning d ON change the aging curve / gamma?

| quantity | value |
| --- | --- |
| curve_corr | 0.999408 |
| curve_max_abs_delta | 0.0109004 |
| peak_age_aging | 2.40862 |
| peak_age_full | 2.52611 |
| gamma_max_abs_delta | 4.41831e-05 |
| gamma_corr | 0.998579 |


### 5. Factor stability across the grid (v_j and u_i)


v_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9311 | 0.9097 | 0.9300 |
| aging | 0.9311 | 1.0000 | 0.8465 | 0.9646 |
| drift | 0.9097 | 0.8465 | 1.0000 | 0.8825 |
| full | 0.9300 | 0.9646 | 0.8825 | 1.0000 |

v_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9104 | 0.8900 | 0.9086 |
| aging | 0.9104 | 1.0000 | 0.7996 | 0.9457 |
| drift | 0.8900 | 0.7996 | 1.0000 | 0.8491 |
| full | 0.9086 | 0.9457 | 0.8491 | 1.0000 |

v_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.02197 | 0.03792 | 0.03418 |
| aging | 0.02197 | 0 | 0.03731 | 0.01949 |
| drift | 0.03792 | 0.03731 | 0 | 0.03622 |
| full | 0.03418 | 0.01949 | 0.03622 | 0 |


u_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9936 | 0.9986 | 0.9923 |
| aging | 0.9936 | 1.0000 | 0.9930 | 0.9989 |
| drift | 0.9986 | 0.9930 | 1.0000 | 0.9937 |
| full | 0.9923 | 0.9989 | 0.9937 | 1.0000 |

u_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9941 | 0.9985 | 0.9927 |
| aging | 0.9941 | 1.0000 | 0.9935 | 0.9989 |
| drift | 0.9985 | 0.9935 | 1.0000 | 0.9942 |
| full | 0.9927 | 0.9989 | 0.9942 | 1.0000 |

u_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 1.012 | 0.2013 | 1.034 |
| aging | 1.012 | 0 | 1.005 | 0.1943 |
| drift | 0.2013 | 1.005 | 0 | 1.026 |
| full | 1.034 | 0.1943 | 1.026 | 0 |


### 6. Aging-curve de-biasing across mrc -- ALL_B (mrc 2 vs 5)

population aging curve, centered to mean 0 on shared A_n in [0, 11.69] yr (100 pts), one APC beta=0 gauge.

peak age (A_n at fastest, yr) per curve:
| mrc | variant | peak_age |
| --- | --- | --- |
| 2 | noD | 2.716 |
| 2 | withD | 2.480 |
| 5 | noD | 2.480 |
| 5 | withD | 2.480 |


curve-distance contrasts (RMS on the shared grid):
| contrast | rms |
| --- | --- |
| gap_noD (everyone vs dedicated, no d) | 0.0075812 |
| gap_withD (everyone vs dedicated, +d) | 0.0041657 |
| d_on_mrc2 (\|\|aging-full\|\| @ everyone) | 0.0074478 |
| d_on_mrc5 (\|\|aging-full\|\| @ dedicated) | 0.0040628 |


(expect gap_withD << gap_noD and d_on_mrc2 >> d_on_mrc5 if no-d/everyone is the contaminated curve.)

## WA_M_14-25_mrc2   nu=8

I=5,371  J=196  N=20,234   cells: baseline / aging / drift / full

### 1. Fit 2x2 (descriptive: loglik / RSS / AIC / BIC)

| cell | label | loglik | rss | aic | bic | n_params |
| --- | --- | --- | --- | --- | --- | --- |
| baseline | no-d / no-aging | 41379.8 | 38.4085 | -71627.6 | -27572.1 | 5567 |
| aging | no-d / +aging | 41884.2 | 36.9631 | -72624.4 | -28521.3 | 5573 |
| drift | +d / no-aging | 45522 | 28.3204 | -75386.9 | -13423 | 8697 |
| full | +d / +aging | 45987.7 | 27.733 | -76317.9 | -14352.2 | 8703 |


loglik on the 2x2 grid:
| loglik | no-aging | +aging |
| --- | --- | --- |
| no-d | 41379.8070 | 41884.1814 |
| +d | 45521.9884 | 45987.7038 |


marginal loglik gains (does a term survive adding the other?):
| effect | d_loglik |
| --- | --- |
| aging_gain \| no-d | 504.3745 |
| aging_gain \| +d | 465.7154 |
| d_gain \| no-aging | 4142.1814 |
| d_gain \| +aging | 4103.5224 |
| interaction (non-additivity) | -38.6590 |


### 2. Contribution overlap in the full fit (aging vs d_i)

| quantity | value |
| --- | --- |
| corr_aging_drift | -0.120368 |
| var_aging | 0.000162454 |
| var_drift | 0.000358425 |
| cov | -2.90468e-05 |
| sd_ratio_d_over_aging | 1.48536 |


(corr near 0 => terms carry distinct signal => both needed; gauge-sensitive in the linear direction.)

### 3. Drift block: does turning aging ON shrink d_i?

| quantity | value |
| --- | --- |
| omega_d2_drift | 0.000125779 |
| omega_d2_full | 0.00011715 |
| n_eligible | 3129 |
| d_sd_drift | 0.0094938 |
| d_sd_full | 0.00911365 |
| frac_improver_drift | 0.520614 |
| frac_improver_full | 0.547779 |
| corr_d_drift_full | 0.934377 |
| max_abs_dd | 0.0417845 |


### 4. Aging block: does turning d ON change the aging curve / gamma?

| quantity | value |
| --- | --- |
| curve_corr | 0.998288 |
| curve_max_abs_delta | 0.0214068 |
| peak_age_aging | 2.69792 |
| peak_age_full | 2.46332 |
| gamma_max_abs_delta | 0.000937164 |
| gamma_corr | 0.988982 |


### 5. Factor stability across the grid (v_j and u_i)


v_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9845 | 0.9333 | 0.9591 |
| aging | 0.9845 | 1.0000 | 0.9288 | 0.9621 |
| drift | 0.9333 | 0.9288 | 1.0000 | 0.9188 |
| full | 0.9591 | 0.9621 | 0.9188 | 1.0000 |

v_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9779 | 0.9231 | 0.9490 |
| aging | 0.9779 | 1.0000 | 0.9214 | 0.9555 |
| drift | 0.9231 | 0.9214 | 1.0000 | 0.9065 |
| full | 0.9490 | 0.9555 | 0.9065 | 1.0000 |

v_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.007921 | 0.02273 | 0.02782 |
| aging | 0.007921 | 0 | 0.0237 | 0.02083 |
| drift | 0.02273 | 0.0237 | 0 | 0.03185 |
| full | 0.02782 | 0.02083 | 0.03185 | 0 |


u_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9897 | 0.9884 | 0.9792 |
| aging | 0.9897 | 1.0000 | 0.9822 | 0.9883 |
| drift | 0.9884 | 0.9822 | 1.0000 | 0.9868 |
| full | 0.9792 | 0.9883 | 0.9868 | 1.0000 |

u_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9886 | 0.9906 | 0.9807 |
| aging | 0.9886 | 1.0000 | 0.9840 | 0.9910 |
| drift | 0.9906 | 0.9840 | 1.0000 | 0.9848 |
| full | 0.9807 | 0.9910 | 0.9848 | 1.0000 |

u_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.2827 | 0.159 | 0.3017 |
| aging | 0.2827 | 0 | 0.2944 | 0.1636 |
| drift | 0.159 | 0.2944 | 0 | 0.3133 |
| full | 0.3017 | 0.1636 | 0.3133 | 0 |


## WA_M_14-25_mrc5   nu=8

I=1,124  J=120  N=8,029   cells: baseline / aging / drift / full

### 1. Fit 2x2 (descriptive: loglik / RSS / AIC / BIC)

| cell | label | loglik | rss | aic | bic | n_params |
| --- | --- | --- | --- | --- | --- | --- |
| baseline | no-d / no-aging | 14008.4 | 25.6556 | -25530.9 | -16841.3 | 1244 |
| aging | no-d / +aging | 14312.8 | 24.2661 | -26127.6 | -17396 | 1250 |
| drift | +d / no-aging | 15849.2 | 18.8479 | -27344.9 | -12128 | 2369 |
| full | +d / +aging | 15915.6 | 18.6617 | -27563.9 | -12647.4 | 2375 |


loglik on the 2x2 grid:
| loglik | no-aging | +aging |
| --- | --- | --- |
| no-d | 14008.4276 | 14312.7786 |
| +d | 15849.1778 | 15915.6495 |


marginal loglik gains (does a term survive adding the other?):
| effect | d_loglik |
| --- | --- |
| aging_gain \| no-d | 304.3509 |
| aging_gain \| +d | 66.4717 |
| d_gain \| no-aging | 1840.7501 |
| d_gain \| +aging | 1602.8709 |
| interaction (non-additivity) | -237.8792 |


### 2. Contribution overlap in the full fit (aging vs d_i)

| quantity | value |
| --- | --- |
| corr_aging_drift | -0.0620933 |
| var_aging | 0.000417541 |
| var_drift | 0.000489294 |
| cov | -2.80694e-05 |
| sd_ratio_d_over_aging | 1.08252 |


(corr near 0 => terms carry distinct signal => both needed; gauge-sensitive in the linear direction.)

### 3. Drift block: does turning aging ON shrink d_i?

| quantity | value |
| --- | --- |
| omega_d2_drift | 0.000142212 |
| omega_d2_full | 9.97483e-05 |
| n_eligible | 1124 |
| d_sd_drift | 0.0108739 |
| d_sd_full | 0.00886475 |
| frac_improver_drift | 0.451957 |
| frac_improver_full | 0.459964 |
| corr_d_drift_full | 0.871218 |
| max_abs_dd | 0.0165986 |


### 4. Aging block: does turning d ON change the aging curve / gamma?

| quantity | value |
| --- | --- |
| curve_corr | 0.999846 |
| curve_max_abs_delta | 0.0119674 |
| peak_age_aging | 2.69349 |
| peak_age_full | 2.75205 |
| gamma_max_abs_delta | 0.000527015 |
| gamma_corr | 0.901766 |


### 5. Factor stability across the grid (v_j and u_i)


v_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9591 | 0.8563 | 0.9239 |
| aging | 0.9591 | 1.0000 | 0.8624 | 0.9787 |
| drift | 0.8563 | 0.8624 | 1.0000 | 0.9097 |
| full | 0.9239 | 0.9787 | 0.9097 | 1.0000 |

v_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9408 | 0.8501 | 0.8986 |
| aging | 0.9408 | 1.0000 | 0.8735 | 0.9745 |
| drift | 0.8501 | 0.8735 | 1.0000 | 0.9227 |
| full | 0.8986 | 0.9745 | 0.9227 | 1.0000 |

v_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.01672 | 0.0429 | 0.02137 |
| aging | 0.01672 | 0 | 0.04911 | 0.01789 |
| drift | 0.0429 | 0.04911 | 0 | 0.03984 |
| full | 0.02137 | 0.01789 | 0.03984 | 0 |


u_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9812 | 0.9872 | 0.9673 |
| aging | 0.9812 | 1.0000 | 0.9787 | 0.9917 |
| drift | 0.9872 | 0.9787 | 1.0000 | 0.9814 |
| full | 0.9673 | 0.9917 | 0.9814 | 1.0000 |

u_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9681 | 0.9861 | 0.9548 |
| aging | 0.9681 | 1.0000 | 0.9714 | 0.9943 |
| drift | 0.9861 | 0.9714 | 1.0000 | 0.9699 |
| full | 0.9548 | 0.9943 | 0.9699 | 1.0000 |

u_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.06274 | 0.14 | 0.1312 |
| aging | 0.06274 | 0 | 0.1565 | 0.1438 |
| drift | 0.14 | 0.1565 | 0 | 0.06139 |
| full | 0.1312 | 0.1438 | 0.06139 | 0 |


### 6. Aging-curve de-biasing across mrc -- WA_M (mrc 2 vs 5)

population aging curve, centered to mean 0 on shared A_n in [0, 11.65] yr (100 pts), one APC beta=0 gauge.

peak age (A_n at fastest, yr) per curve:
| mrc | variant | peak_age |
| --- | --- | --- |
| 2 | noD | 2.707 |
| 2 | withD | 2.472 |
| 5 | noD | 2.707 |
| 5 | withD | 2.825 |


curve-distance contrasts (RMS on the shared grid):
| contrast | rms |
| --- | --- |
| gap_noD (everyone vs dedicated, no d) | 0.0056553 |
| gap_withD (everyone vs dedicated, +d) | 0.0031348 |
| d_on_mrc2 (\|\|aging-full\|\| @ everyone) | 0.0078262 |
| d_on_mrc5 (\|\|aging-full\|\| @ dedicated) | 0.0047382 |


(expect gap_withD << gap_noD and d_on_mrc2 >> d_on_mrc5 if no-d/everyone is the contaminated curve.)

## WA_W_14-25_mrc2   nu=8

I=3,431  J=150  N=12,751   cells: baseline / aging / drift / full

### 1. Fit 2x2 (descriptive: loglik / RSS / AIC / BIC)

| cell | label | loglik | rss | aic | bic | n_params |
| --- | --- | --- | --- | --- | --- | --- |
| baseline | no-d / no-aging | 23774.4 | 29.0914 | -40388.9 | -13705.8 | 3581 |
| aging | no-d / +aging | 24231.4 | 27.2867 | -41290.7 | -14563 | 3587 |
| drift | +d / no-aging | 26847 | 19.4587 | -43576.8 | -5872.97 | 5575 |
| full | +d / +aging | 27122.3 | 18.9653 | -44151.2 | -6536.6 | 5581 |


loglik on the 2x2 grid:
| loglik | no-aging | +aging |
| --- | --- | --- |
| no-d | 23774.4289 | 24231.3660 |
| +d | 26847.0243 | 27122.2705 |


marginal loglik gains (does a term survive adding the other?):
| effect | d_loglik |
| --- | --- |
| aging_gain \| no-d | 456.9371 |
| aging_gain \| +d | 275.2462 |
| d_gain \| no-aging | 3072.5954 |
| d_gain \| +aging | 2890.9044 |
| interaction (non-additivity) | -181.6909 |


### 2. Contribution overlap in the full fit (aging vs d_i)

| quantity | value |
| --- | --- |
| corr_aging_drift | -0.0727792 |
| var_aging | 0.000330507 |
| var_drift | 0.000503809 |
| cov | -2.97005e-05 |
| sd_ratio_d_over_aging | 1.23465 |


(corr near 0 => terms carry distinct signal => both needed; gauge-sensitive in the linear direction.)

### 3. Drift block: does turning aging ON shrink d_i?

| quantity | value |
| --- | --- |
| omega_d2_drift | 0.000175604 |
| omega_d2_full | 0.000156491 |
| n_eligible | 1993 |
| d_sd_drift | 0.0114073 |
| d_sd_full | 0.0106721 |
| frac_improver_drift | 0.487205 |
| frac_improver_full | 0.522328 |
| corr_d_drift_full | 0.921834 |
| max_abs_dd | 0.0282087 |


### 4. Aging block: does turning d ON change the aging curve / gamma?

| quantity | value |
| --- | --- |
| curve_corr | 0.999101 |
| curve_max_abs_delta | 0.0197721 |
| peak_age_aging | 2.69792 |
| peak_age_full | 2.63927 |
| gamma_max_abs_delta | 0.000870861 |
| gamma_corr | 0.314911 |


### 5. Factor stability across the grid (v_j and u_i)


v_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9848 | 0.9510 | 0.9775 |
| aging | 0.9848 | 1.0000 | 0.9493 | 0.9842 |
| drift | 0.9510 | 0.9493 | 1.0000 | 0.9614 |
| full | 0.9775 | 0.9842 | 0.9614 | 1.0000 |

v_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9802 | 0.9595 | 0.9758 |
| aging | 0.9802 | 1.0000 | 0.9569 | 0.9823 |
| drift | 0.9595 | 0.9569 | 1.0000 | 0.9719 |
| full | 0.9758 | 0.9823 | 0.9719 | 1.0000 |

v_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.01406 | 0.03031 | 0.01558 |
| aging | 0.01406 | 0 | 0.02762 | 0.01762 |
| drift | 0.03031 | 0.02762 | 0 | 0.03158 |
| full | 0.01558 | 0.01762 | 0.03158 | 0 |


u_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9883 | 0.9880 | 0.9821 |
| aging | 0.9883 | 1.0000 | 0.9840 | 0.9945 |
| drift | 0.9880 | 0.9840 | 1.0000 | 0.9839 |
| full | 0.9821 | 0.9945 | 0.9839 | 1.0000 |

u_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9835 | 0.9881 | 0.9751 |
| aging | 0.9835 | 1.0000 | 0.9777 | 0.9932 |
| drift | 0.9881 | 0.9777 | 1.0000 | 0.9753 |
| full | 0.9751 | 0.9932 | 0.9753 | 1.0000 |

u_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.2075 | 0.2115 | 0.2215 |
| aging | 0.2075 | 0 | 0.1824 | 0.1248 |
| drift | 0.2115 | 0.1824 | 0 | 0.1843 |
| full | 0.2215 | 0.1248 | 0.1843 | 0 |


## WA_W_14-25_mrc5   nu=8

I=646  J=69  N=4,441   cells: baseline / aging / drift / full

### 1. Fit 2x2 (descriptive: loglik / RSS / AIC / BIC)

| cell | label | loglik | rss | aic | bic | n_params |
| --- | --- | --- | --- | --- | --- | --- |
| baseline | no-d / no-aging | 6854.51 | 17.2252 | -12281 | -7712.39 | 715 |
| aging | no-d / +aging | 7071.52 | 15.895 | -12703 | -8096.02 | 721 |
| drift | +d / no-aging | 8046.63 | 11.2322 | -13585.5 | -5562.35 | 1362 |
| full | +d / +aging | 8050.16 | 11.2101 | -13635 | -5747.55 | 1368 |


loglik on the 2x2 grid:
| loglik | no-aging | +aging |
| --- | --- | --- |
| no-d | 6854.5095 | 7071.5163 |
| +d | 8046.6281 | 8050.1551 |


marginal loglik gains (does a term survive adding the other?):
| effect | d_loglik |
| --- | --- |
| aging_gain \| no-d | 217.0068 |
| aging_gain \| +d | 3.5270 |
| d_gain \| no-aging | 1192.1186 |
| d_gain \| +aging | 978.6389 |
| interaction (non-additivity) | -213.4797 |


### 2. Contribution overlap in the full fit (aging vs d_i)

| quantity | value |
| --- | --- |
| corr_aging_drift | 0.00285885 |
| var_aging | 0.00062989 |
| var_drift | 0.000746162 |
| cov | 1.96037e-06 |
| sd_ratio_d_over_aging | 1.08839 |


(corr near 0 => terms carry distinct signal => both needed; gauge-sensitive in the linear direction.)

### 3. Drift block: does turning aging ON shrink d_i?

| quantity | value |
| --- | --- |
| omega_d2_drift | 0.000201847 |
| omega_d2_full | 0.000145066 |
| n_eligible | 646 |
| d_sd_drift | 0.0129981 |
| d_sd_full | 0.010738 |
| frac_improver_drift | 0.47678 |
| frac_improver_full | 0.467492 |
| corr_d_drift_full | 0.896311 |
| max_abs_dd | 0.0222921 |


### 4. Aging block: does turning d ON change the aging curve / gamma?

| quantity | value |
| --- | --- |
| curve_corr | 0.998004 |
| curve_max_abs_delta | 0.0040184 |
| peak_age_aging | 2.75205 |
| peak_age_full | 3.16193 |
| gamma_max_abs_delta | 0.000177107 |
| gamma_corr | 0.913508 |


### 5. Factor stability across the grid (v_j and u_i)


v_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9531 | 0.9077 | 0.9540 |
| aging | 0.9531 | 1.0000 | 0.9039 | 0.9897 |
| drift | 0.9077 | 0.9039 | 1.0000 | 0.9380 |
| full | 0.9540 | 0.9897 | 0.9380 | 1.0000 |

v_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9332 | 0.9261 | 0.9360 |
| aging | 0.9332 | 1.0000 | 0.9269 | 0.9875 |
| drift | 0.9261 | 0.9269 | 1.0000 | 0.9556 |
| full | 0.9360 | 0.9875 | 0.9556 | 1.0000 |

v_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.02451 | 0.04324 | 0.0234 |
| aging | 0.02451 | 0 | 0.04212 | 0.01488 |
| drift | 0.04324 | 0.04212 | 0 | 0.03444 |
| full | 0.0234 | 0.01488 | 0.03444 | 0 |


u_j Pearson corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9723 | 0.9922 | 0.9719 |
| aging | 0.9723 | 1.0000 | 0.9730 | 0.9964 |
| drift | 0.9922 | 0.9730 | 1.0000 | 0.9789 |
| full | 0.9719 | 0.9964 | 0.9789 | 1.0000 |

u_j Spearman corr (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 1.0000 | 0.9494 | 0.9880 | 0.9504 |
| aging | 0.9494 | 1.0000 | 0.9526 | 0.9950 |
| drift | 0.9880 | 0.9526 | 1.0000 | 0.9621 |
| full | 0.9504 | 0.9950 | 0.9621 | 1.0000 |

u_j max|delta| (4x4):
| cell | baseline | aging | drift | full |
| --- | --- | --- | --- | --- |
| baseline | 0 | 0.09151 | 0.06913 | 0.09906 |
| aging | 0.09151 | 0 | 0.1119 | 0.07484 |
| drift | 0.06913 | 0.1119 | 0 | 0.08485 |
| full | 0.09906 | 0.07484 | 0.08485 | 0 |


### 6. Aging-curve de-biasing across mrc -- WA_W (mrc 2 vs 5)

population aging curve, centered to mean 0 on shared A_n in [0, 11.65] yr (100 pts), one APC beta=0 gauge.

peak age (A_n at fastest, yr) per curve:
| mrc | variant | peak_age |
| --- | --- | --- |
| 2 | noD | 2.707 |
| 2 | withD | 2.589 |
| 5 | noD | 2.825 |
| 5 | withD | 3.178 |


curve-distance contrasts (RMS on the shared grid):
| contrast | rms |
| --- | --- |
| gap_noD (everyone vs dedicated, no d) | 0.0065148 |
| gap_withD (everyone vs dedicated, +d) | 0.0024907 |
| d_on_mrc2 (\|\|aging-full\|\| @ everyone) | 0.0070616 |
| d_on_mrc5 (\|\|aging-full\|\| @ dedicated) | 0.0018669 |


(expect gap_withD << gap_noD and d_on_mrc2 >> d_on_mrc5 if no-d/everyone is the contaminated curve.)

CSV rollups -> results/model_selection/aging_vs_drift
  fit_grid.csv  gains.csv  block_shift.csv  stability.csv  curve_debias.csv
