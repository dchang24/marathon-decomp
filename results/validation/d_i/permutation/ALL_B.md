# d_i v_j-bias permutation null -- ALL_B

d_i v_j-bias permutation null -- ALL_B   nu=8   n_perm=300
H0: within-athlete career-stage unrelated to performance. Expect the no-d statistic in the upper tail, the +d statistic in the null.

## ALL_B_14-25_mrc2   nu=8   N=1,244,290  J=347   (+300 draws; had 200)

| fit | n_perm_total | T_obs | null_mean | null_sd | z | p_one | p_two |
| --- | --- | --- | --- | --- | --- | --- | --- |
| no-d | 500 | 0.0577 | -0.0893 | 0.0416 | 3.5372 | 0.0020 | 0.7824 |
| +d | 500 | -0.0233 | -0.0593 | 0.0409 | 0.8792 | 0.1936 | 0.8303 |


no-d significant (p_one=0.001996, n=500), +d null as expected

## ALL_B_14-25_mrc5   nu=8   N=432,253  J=330   (+300 draws; had 200)

| fit | n_perm_total | T_obs | null_mean | null_sd | z | p_one | p_two |
| --- | --- | --- | --- | --- | --- | --- | --- |
| no-d | 500 | 0.2318 | -0.0278 | 0.0540 | 4.8072 | 0.0020 | 0.0020 |
| +d | 500 | 0.0726 | -0.0187 | 0.0526 | 1.7347 | 0.0579 | 0.2136 |


no-d significant (p_one=0.001996, n=500), +d null as expected

## ALL_B summary (p over full accumulated pool)

| slug | fit | n_perm_total | T_obs | z | p_one | p_two |
| --- | --- | --- | --- | --- | --- | --- |
| ALL_B_14-25_mrc2 | +d | 500 | -0.0233 | 0.8792 | 0.1936 | 0.8303 |
| ALL_B_14-25_mrc2 | no-d | 500 | 0.0577 | 3.5372 | 0.0020 | 0.7824 |
| ALL_B_14-25_mrc5 | +d | 500 | 0.0726 | 1.7347 | 0.0579 | 0.2136 |
| ALL_B_14-25_mrc5 | no-d | 500 | 0.2318 | 4.8072 | 0.0020 | 0.0020 |

