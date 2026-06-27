# d_i v_j-bias permutation null -- ALL_M

d_i v_j-bias permutation null -- ALL_M   nu=8   n_perm=300
H0: within-athlete career-stage unrelated to performance. Expect the no-d statistic in the upper tail, the +d statistic in the null.

## ALL_M_14-25_mrc2   nu=8   N=863,072  J=344   (+300 draws; had 200)

| fit | n_perm_total | T_obs | null_mean | null_sd | z | p_one | p_two |
| --- | --- | --- | --- | --- | --- | --- | --- |
| no-d | 500 | 0.0557 | -0.0970 | 0.0422 | 3.6160 | 0.0020 | 0.8363 |
| +d | 500 | -0.0290 | -0.0716 | 0.0397 | 1.0760 | 0.1497 | 0.8523 |


no-d significant (p_one=0.001996, n=500), +d null as expected

## ALL_M_14-25_mrc5   nu=8   N=299,892  J=328   (+300 draws; had 200)

| fit | n_perm_total | T_obs | null_mean | null_sd | z | p_one | p_two |
| --- | --- | --- | --- | --- | --- | --- | --- |
| no-d | 500 | 0.2054 | -0.0284 | 0.0559 | 4.1818 | 0.0020 | 0.0020 |
| +d | 500 | 0.0374 | -0.0223 | 0.0547 | 1.0900 | 0.1437 | 0.5289 |


no-d significant (p_one=0.001996, n=500), +d null as expected

## ALL_M summary (p over full accumulated pool)

| slug | fit | n_perm_total | T_obs | z | p_one | p_two |
| --- | --- | --- | --- | --- | --- | --- |
| ALL_M_14-25_mrc2 | +d | 500 | -0.0290 | 1.0760 | 0.1497 | 0.8523 |
| ALL_M_14-25_mrc2 | no-d | 500 | 0.0557 | 3.6160 | 0.0020 | 0.8363 |
| ALL_M_14-25_mrc5 | +d | 500 | 0.0374 | 1.0900 | 0.1437 | 0.5289 |
| ALL_M_14-25_mrc5 | no-d | 500 | 0.2054 | 4.1818 | 0.0020 | 0.0020 |

