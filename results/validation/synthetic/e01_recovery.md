
# e01 -- full-model recovery (Test A)


## Per-seed rows

160 fits over designs=['base', 'tiny'], samplers=['staggered', 'balanced'], nus=[inf, 6.0], seeds=20.

## Recovery (mean across seeds; r_v/r_u gauge-fixed, raw shown too)

| design | sampler | nu | r_v__mean | r_v_raw__mean | r_u__mean | r_d__mean | sigma2_ratio__mean | omega_d2_ratio__mean | mre__mean |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| base | balanced | 6 | 0.9931 | 0.9896 | 0.9993 | - | 0.944 | 0.05208 | 0.0268 |
| base | balanced | - | 0.9944 | 0.9911 | 0.9994 | - | 0.9414 | 0.05208 | 0.02322 |
| base | staggered | 6 | 0.9786 | 0.9769 | 0.9993 | - | 0.9321 | 0.05208 | 0.02664 |
| base | staggered | - | 0.9867 | 0.9861 | 0.9994 | - | 0.9472 | 0.05208 | 0.02332 |
| tiny | balanced | 6 | 0.9786 | 0.9575 | 0.999 | - | 0.9162 | 0.05208 | 0.02642 |
| tiny | balanced | - | 0.9845 | 0.972 | 0.9993 | - | 0.9248 | 0.05208 | 0.02303 |
| tiny | staggered | 6 | 0.9624 | 0.9565 | 0.9989 | - | 0.9235 | 0.05208 | 0.02651 |
| tiny | staggered | - | 0.9811 | 0.9781 | 0.9992 | - | 0.9332 | 0.05208 | 0.02312 |


## Aging curvature error (full model confounds with drift; no-drift refit is the identified value; true=0.0044)

| design | sampler | nu | curvature_err__mean | curvature_err_nodrift__mean |
| --- | --- | --- | --- | --- |
| base | balanced | 6 | -15.93 | -15.93 |
| base | balanced | - | 31.23 | 31.23 |
| base | staggered | 6 | -10.38 | -10.38 |
| base | staggered | - | 28.67 | 28.67 |
| tiny | balanced | 6 | 281.7 | 281.7 |
| tiny | balanced | - | -287.3 | -287.3 |
| tiny | staggered | 6 | -20.28 | -20.28 |
| tiny | staggered | - | 98.15 | 98.15 |


## Convergence (mean across seeds)

| design | sampler | nu | converged__mean | n_iter__mean | stationarity_v__mean | oracle_margin__mean | mono_min_step__mean |
| --- | --- | --- | --- | --- | --- | --- | --- |
| base | balanced | 6 | 1 | 10.45 | 3.814e-06 | 223.5 | 0.0005248 |
| base | balanced | - | 1 | 7.4 | 2.274e-06 | 227 | 0.0002248 |
| base | staggered | 6 | 1 | 7 | 5.029e-06 | 225.6 | 0.0003027 |
| base | staggered | - | 1 | 5.5 | 2.525e-06 | 227.8 | 0.0003747 |
| tiny | balanced | 6 | 1 | 11.15 | 1.038e-06 | 91.72 | 0.0001318 |
| tiny | balanced | - | 1 | 8.35 | 5.329e-07 | 94.97 | 9.67e-05 |
| tiny | staggered | 6 | 1 | 7.1 | 4.493e-06 | 95.28 | 0.0001194 |
| tiny | staggered | - | 1 | 5.45 | 5.083e-06 | 91.97 | 0.0001344 |


## ALS == Anderson at the fixed point

| design | nu | seed | als_n_iter | and_n_iter | agree_r_v | sigma2_ratio |
| --- | --- | --- | --- | --- | --- | --- |
| base | - | 2025 | 55 | 8 | 0.9996 | 1 |
| tiny | - | 2025 | 70 | 8 | 0.9997 | 1 |

