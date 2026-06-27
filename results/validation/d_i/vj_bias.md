# d_i v_j-bias test

d_i v_j-bias test   nu=8   mrc=[2, 5]   slices: Po10_W, Po10_M, Po10_B, ALL_W, ALL_M, ALL_B
delta tables present: ['eb', 'loo']
test: corr(delta, v_no-d) > 0  and  corr(delta, v_+d) ~ 0 (year-partialled). Significance deferred to q04 permutation null.

## Po10_W_14-25_mrc2   nu=8

J=203 races   v: no-d=agingS4gv, +d=full   estimators present: ['eb', 'loo']

delta agreement  corr(delta_EB, delta_LOO) = 0.4334

### estimator = EB

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.3105 | -0.0682 |
| spearman_raw | 0.3512 | -0.0714 |
| pearson_partial | 0.2258 | 0.0139 |
| spearman_partial | 0.3331 | 0.1057 |
| slope_partial | 1.1041 | 0.0657 |


year-partialled pearson:  no-d=+0.2258  +d=+0.0139  drop=+0.2119  frac_removed=0.938  -> PASS

### estimator = LOO

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.1900 | 0.0479 |
| spearman_raw | 0.2698 | 0.1408 |
| pearson_partial | 0.2099 | 0.0416 |
| spearman_partial | 0.3056 | 0.1342 |
| slope_partial | 0.8300 | 0.1593 |


year-partialled pearson:  no-d=+0.2099  +d=+0.0416  drop=+0.1682  frac_removed=0.802  -> PASS

## Po10_W_14-25_mrc5   nu=8

J=166 races   v: no-d=agingS4gv, +d=full   estimators present: ['eb', 'loo']

delta agreement  corr(delta_EB, delta_LOO) = 0.7325

### estimator = EB

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.3439 | 0.0832 |
| spearman_raw | 0.4259 | 0.1548 |
| pearson_partial | 0.2366 | -0.0221 |
| spearman_partial | 0.3417 | 0.0747 |
| slope_partial | 1.0068 | -0.0909 |


year-partialled pearson:  no-d=+0.2366  +d=-0.0221  drop=+0.2587  frac_removed=1.094  -> PASS

### estimator = LOO

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.1289 | -0.0650 |
| spearman_raw | 0.2062 | 0.0163 |
| pearson_partial | 0.2203 | -0.0130 |
| spearman_partial | 0.3229 | 0.0886 |
| slope_partial | 0.7433 | -0.0422 |


year-partialled pearson:  no-d=+0.2203  +d=-0.0130  drop=+0.2333  frac_removed=1.059  -> PASS

## Po10_M_14-25_mrc2   nu=8

J=252 races   v: no-d=agingS4gv, +d=full   estimators present: ['eb', 'loo']

delta agreement  corr(delta_EB, delta_LOO) = 0.4747

### estimator = EB

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.4132 | 0.0878 |
| spearman_raw | 0.4941 | 0.1666 |
| pearson_partial | 0.2271 | 0.0563 |
| spearman_partial | 0.3272 | 0.1546 |
| slope_partial | 1.2789 | 0.3022 |


year-partialled pearson:  no-d=+0.2271  +d=+0.0563  drop=+0.1709  frac_removed=0.752  -> PASS

### estimator = LOO

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.3469 | 0.2177 |
| spearman_raw | 0.3743 | 0.2873 |
| pearson_partial | 0.3410 | 0.2133 |
| spearman_partial | 0.4044 | 0.2841 |
| slope_partial | 1.3515 | 0.8064 |


year-partialled pearson:  no-d=+0.3410  +d=+0.2133  drop=+0.1276  frac_removed=0.374  -> check

## Po10_M_14-25_mrc5   nu=8

J=216 races   v: no-d=agingS4gv, +d=full   estimators present: ['eb', 'loo']

delta agreement  corr(delta_EB, delta_LOO) = 0.7921

### estimator = EB

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.4972 | 0.2904 |
| spearman_raw | 0.5759 | 0.4029 |
| pearson_partial | 0.2720 | 0.0476 |
| spearman_partial | 0.3695 | 0.1649 |
| slope_partial | 1.1493 | 0.1917 |


year-partialled pearson:  no-d=+0.2720  +d=+0.0476  drop=+0.2244  frac_removed=0.825  -> PASS

### estimator = LOO

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.2893 | 0.1217 |
| spearman_raw | 0.3173 | 0.1753 |
| pearson_partial | 0.3227 | 0.1157 |
| spearman_partial | 0.3920 | 0.2044 |
| slope_partial | 1.0751 | 0.3674 |


year-partialled pearson:  no-d=+0.3227  +d=+0.1157  drop=+0.2070  frac_removed=0.642  -> PASS

## Po10_B_14-25_mrc2   nu=8

J=265 races   v: no-d=agingS4gv, +d=full   estimators present: ['eb', 'loo']

delta agreement  corr(delta_EB, delta_LOO) = 0.4322

### estimator = EB

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.3447 | -0.0132 |
| spearman_raw | 0.4024 | 0.0142 |
| pearson_partial | 0.2105 | 0.0330 |
| spearman_partial | 0.2671 | 0.0720 |
| slope_partial | 1.2234 | 0.1847 |


year-partialled pearson:  no-d=+0.2105  +d=+0.0330  drop=+0.1775  frac_removed=0.843  -> PASS

### estimator = LOO

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.3004 | 0.1652 |
| spearman_raw | 0.3598 | 0.2388 |
| pearson_partial | 0.3001 | 0.1671 |
| spearman_partial | 0.3871 | 0.2406 |
| slope_partial | 1.3323 | 0.7140 |


year-partialled pearson:  no-d=+0.3001  +d=+0.1671  drop=+0.1330  frac_removed=0.443  -> check

## Po10_B_14-25_mrc5   nu=8

J=236 races   v: no-d=agingS4gv, +d=full   estimators present: ['eb', 'loo']

delta agreement  corr(delta_EB, delta_LOO) = 0.7686

### estimator = EB

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.4951 | 0.2580 |
| spearman_raw | 0.5384 | 0.3472 |
| pearson_partial | 0.2882 | 0.0457 |
| spearman_partial | 0.3455 | 0.1388 |
| slope_partial | 1.1908 | 0.1795 |


year-partialled pearson:  no-d=+0.2882  +d=+0.0457  drop=+0.2425  frac_removed=0.841  -> PASS

### estimator = LOO

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.3012 | 0.1189 |
| spearman_raw | 0.2814 | 0.1458 |
| pearson_partial | 0.3257 | 0.1104 |
| spearman_partial | 0.3583 | 0.1717 |
| slope_partial | 1.1031 | 0.3554 |


year-partialled pearson:  no-d=+0.3257  +d=+0.1104  drop=+0.2153  frac_removed=0.661  -> PASS

## ALL_W_14-25_mrc2   nu=8

J=337 races   v: no-d=agingS4gv, +d=full   estimators present: ['eb', 'loo']

delta agreement  corr(delta_EB, delta_LOO) = 0.4291

### estimator = EB

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.0560 | -0.1680 |
| spearman_raw | 0.0612 | -0.1690 |
| pearson_partial | 0.1860 | 0.0351 |
| spearman_partial | 0.2169 | 0.0704 |
| slope_partial | 1.2944 | 0.2390 |


year-partialled pearson:  no-d=+0.1860  +d=+0.0351  drop=+0.1509  frac_removed=0.811  -> PASS

### estimator = LOO

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.0537 | -0.0225 |
| spearman_raw | 0.0435 | -0.0297 |
| pearson_partial | 0.0511 | -0.0308 |
| spearman_partial | 0.0409 | -0.0317 |
| slope_partial | 0.2957 | -0.1744 |


year-partialled pearson:  no-d=+0.0511  +d=-0.0308  drop=+0.0820  frac_removed=1.603  -> check

## ALL_W_14-25_mrc5   nu=8

J=295 races   v: no-d=agingS4gv, +d=full   estimators present: ['eb', 'loo']

delta agreement  corr(delta_EB, delta_LOO) = 0.6685

### estimator = EB

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.2421 | -0.0226 |
| spearman_raw | 0.2625 | -0.0173 |
| pearson_partial | 0.3113 | 0.0952 |
| spearman_partial | 0.3499 | 0.1520 |
| slope_partial | 1.4983 | 0.4314 |


year-partialled pearson:  no-d=+0.3113  +d=+0.0952  drop=+0.2161  frac_removed=0.694  -> PASS

### estimator = LOO

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.2531 | 0.0793 |
| spearman_raw | 0.2736 | 0.0990 |
| pearson_partial | 0.2500 | 0.0623 |
| spearman_partial | 0.2707 | 0.0987 |
| slope_partial | 0.9551 | 0.2240 |


year-partialled pearson:  no-d=+0.2500  +d=+0.0623  drop=+0.1877  frac_removed=0.751  -> PASS

## ALL_M_14-25_mrc2   nu=8

J=344 races   v: no-d=agingS4gv, +d=full   estimators present: ['eb', 'loo']

delta agreement  corr(delta_EB, delta_LOO) = 0.3999

### estimator = EB

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.1040 | -0.1119 |
| spearman_raw | 0.0652 | -0.1681 |
| pearson_partial | 0.2202 | 0.0858 |
| spearman_partial | 0.2490 | 0.1180 |
| slope_partial | 1.6222 | 0.6122 |


year-partialled pearson:  no-d=+0.2202  +d=+0.0858  drop=+0.1344  frac_removed=0.610  -> PASS

### estimator = LOO

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.0631 | -0.0063 |
| spearman_raw | 0.0696 | -0.0044 |
| pearson_partial | 0.0557 | -0.0290 |
| spearman_partial | 0.0665 | -0.0092 |
| slope_partial | 0.3548 | -0.1786 |


year-partialled pearson:  no-d=+0.0557  +d=-0.0290  drop=+0.0847  frac_removed=1.520  -> check

## ALL_M_14-25_mrc5   nu=8

J=328 races   v: no-d=agingS4gv, +d=full   estimators present: ['eb', 'loo']

delta agreement  corr(delta_EB, delta_LOO) = 0.7620

### estimator = EB

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.2016 | -0.0386 |
| spearman_raw | 0.2247 | -0.0095 |
| pearson_partial | 0.2289 | 0.0406 |
| spearman_partial | 0.2681 | 0.0963 |
| slope_partial | 1.2538 | 0.2133 |


year-partialled pearson:  no-d=+0.2289  +d=+0.0406  drop=+0.1883  frac_removed=0.823  -> PASS

### estimator = LOO

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.2055 | 0.0358 |
| spearman_raw | 0.2568 | 0.1029 |
| pearson_partial | 0.2054 | 0.0374 |
| spearman_partial | 0.2582 | 0.1045 |
| slope_partial | 0.9921 | 0.1731 |


year-partialled pearson:  no-d=+0.2054  +d=+0.0374  drop=+0.1681  frac_removed=0.818  -> PASS

## ALL_B_14-25_mrc2   nu=8

J=347 races   v: no-d=agingS4gv, +d=full   estimators present: ['eb', 'loo']

delta agreement  corr(delta_EB, delta_LOO) = 0.3962

### estimator = EB

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.0747 | -0.1453 |
| spearman_raw | 0.0332 | -0.1954 |
| pearson_partial | 0.2038 | 0.0643 |
| spearman_partial | 0.2189 | 0.0858 |
| slope_partial | 1.5023 | 0.4611 |


year-partialled pearson:  no-d=+0.2038  +d=+0.0643  drop=+0.1395  frac_removed=0.685  -> PASS

### estimator = LOO

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.0657 | -0.0012 |
| spearman_raw | 0.0382 | -0.0158 |
| pearson_partial | 0.0577 | -0.0233 |
| spearman_partial | 0.0326 | -0.0234 |
| slope_partial | 0.3664 | -0.1442 |


year-partialled pearson:  no-d=+0.0577  +d=-0.0233  drop=+0.0810  frac_removed=1.405  -> PASS

## ALL_B_14-25_mrc5   nu=8

J=330 races   v: no-d=agingS4gv, +d=full   estimators present: ['eb', 'loo']

delta agreement  corr(delta_EB, delta_LOO) = 0.7243

### estimator = EB

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.2045 | -0.0322 |
| spearman_raw | 0.2244 | -0.0182 |
| pearson_partial | 0.2550 | 0.0738 |
| spearman_partial | 0.2966 | 0.1286 |
| slope_partial | 1.4489 | 0.4001 |


year-partialled pearson:  no-d=+0.2550  +d=+0.0738  drop=+0.1813  frac_removed=0.711  -> PASS

### estimator = LOO

| quantity | no-d (agingS4gv) | +d (full) |
| --- | --- | --- |
| pearson_raw | 0.2319 | 0.0740 |
| spearman_raw | 0.2637 | 0.1137 |
| pearson_partial | 0.2318 | 0.0726 |
| spearman_partial | 0.2634 | 0.1173 |
| slope_partial | 1.1868 | 0.3547 |


year-partialled pearson:  no-d=+0.2318  +d=+0.0726  drop=+0.1592  frac_removed=0.687  -> PASS

## Headline: year-partialled pearson across slices

| slug | estimator | corr_nod | corr_d | drop | frac_removed | pass_heuristic |
| --- | --- | --- | --- | --- | --- | --- |
| Po10_W_14-25_mrc2 | eb | 0.2258 | 0.0139 | 0.2119 | 0.9385 | True |
| Po10_W_14-25_mrc2 | loo | 0.2099 | 0.0416 | 0.1682 | 0.8016 | True |
| Po10_W_14-25_mrc5 | eb | 0.2366 | -0.0221 | 0.2587 | 1.0936 | True |
| Po10_W_14-25_mrc5 | loo | 0.2203 | -0.0130 | 0.2333 | 1.0588 | True |
| Po10_M_14-25_mrc2 | eb | 0.2271 | 0.0563 | 0.1709 | 0.7523 | True |
| Po10_M_14-25_mrc2 | loo | 0.3410 | 0.2133 | 0.1276 | 0.3743 | False |
| Po10_M_14-25_mrc5 | eb | 0.2720 | 0.0476 | 0.2244 | 0.8250 | True |
| Po10_M_14-25_mrc5 | loo | 0.3227 | 0.1157 | 0.2070 | 0.6415 | True |
| Po10_B_14-25_mrc2 | eb | 0.2105 | 0.0330 | 0.1775 | 0.8431 | True |
| Po10_B_14-25_mrc2 | loo | 0.3001 | 0.1671 | 0.1330 | 0.4431 | False |
| Po10_B_14-25_mrc5 | eb | 0.2882 | 0.0457 | 0.2425 | 0.8415 | True |
| Po10_B_14-25_mrc5 | loo | 0.3257 | 0.1104 | 0.2153 | 0.6611 | True |
| ALL_W_14-25_mrc2 | eb | 0.1860 | 0.0351 | 0.1509 | 0.8112 | True |
| ALL_W_14-25_mrc2 | loo | 0.0511 | -0.0308 | 0.0820 | 1.6032 | False |
| ALL_W_14-25_mrc5 | eb | 0.3113 | 0.0952 | 0.2161 | 0.6942 | True |
| ALL_W_14-25_mrc5 | loo | 0.2500 | 0.0623 | 0.1877 | 0.7509 | True |
| ALL_M_14-25_mrc2 | eb | 0.2202 | 0.0858 | 0.1344 | 0.6104 | True |
| ALL_M_14-25_mrc2 | loo | 0.0557 | -0.0290 | 0.0847 | 1.5198 | False |
| ALL_M_14-25_mrc5 | eb | 0.2289 | 0.0406 | 0.1883 | 0.8227 | True |
| ALL_M_14-25_mrc5 | loo | 0.2054 | 0.0374 | 0.1681 | 0.8182 | True |
| ALL_B_14-25_mrc2 | eb | 0.2038 | 0.0643 | 0.1395 | 0.6845 | True |
| ALL_B_14-25_mrc2 | loo | 0.0577 | -0.0233 | 0.0810 | 1.4045 | True |
| ALL_B_14-25_mrc5 | eb | 0.2550 | 0.0738 | 0.1813 | 0.7108 | True |
| ALL_B_14-25_mrc5 | loo | 0.2318 | 0.0726 | 0.1592 | 0.6869 | True |


CSV -> results/validation/d_i
  bias_corr.csv  bias_summary.csv  delta_agreement.csv
