# d_i v_j-bias permutation null -- ALL_W

d_i v_j-bias permutation null -- ALL_W   nu=8   n_perm=300
H0: within-athlete career-stage unrelated to performance. Expect the no-d statistic in the upper tail, the +d statistic in the null.

## ALL_W_14-25_mrc2   nu=8   N=381,031  J=337   (+300 draws; had 200)

| fit | n_perm_total | T_obs | null_mean | null_sd | z | p_one | p_two |
| --- | --- | --- | --- | --- | --- | --- | --- |
| no-d | 500 | 0.0511 | -0.0476 | 0.0488 | 2.0225 | 0.0200 | 0.4990 |
| +d | 500 | -0.0308 | -0.0323 | 0.0521 | 0.0283 | 0.4810 | 0.6427 |


no-d significant (p_one=0.01996, n=500), +d null as expected

## ALL_W_14-25_mrc5   nu=8   N=131,398  J=295   (+300 draws; had 200)

| fit | n_perm_total | T_obs | null_mean | null_sd | z | p_one | p_two |
| --- | --- | --- | --- | --- | --- | --- | --- |
| no-d | 500 | 0.2500 | -0.0188 | 0.0633 | 4.2432 | 0.0020 | 0.0020 |
| +d | 500 | 0.0623 | -0.0242 | 0.0628 | 1.3770 | 0.0878 | 0.3713 |


no-d significant (p_one=0.001996, n=500), +d null as expected

## ALL_W summary (p over full accumulated pool)

| slug | fit | n_perm_total | T_obs | z | p_one | p_two |
| --- | --- | --- | --- | --- | --- | --- |
| ALL_W_14-25_mrc2 | +d | 500 | -0.0308 | 0.0283 | 0.4810 | 0.6427 |
| ALL_W_14-25_mrc2 | no-d | 500 | 0.0511 | 2.0225 | 0.0200 | 0.4990 |
| ALL_W_14-25_mrc5 | +d | 500 | 0.0623 | 1.3770 | 0.0878 | 0.3713 |
| ALL_W_14-25_mrc5 | no-d | 500 | 0.2500 | 4.2432 | 0.0020 | 0.0020 |

