# Identity-mismatch sensitivity -- Po10_B_14-25_mrc2/full_nu8p00_best__7cde3824

Inject ADDITIONAL linkage error on top of the production resolution -- `op` = `break` (split one runner's records into two, mimicking a failed join), `join` (merge two different runners, mimicking a wrong match), or `both` -- at requested rate `p`, refit under Monte-Carlo replication, and compare the perturbed race factor `v_j` to the unperturbed baseline on the common race set (rank 1 = hardest).

**sens_ratio** = mismatch per-race SD / bootstrap per-race SD: the mismatch-induced movement relative to the athlete-sampling noise the analysis already carries. `sens_ratio < 1` means identity error perturbs `v_j` by LESS than ordinary sampling noise. **p\*** = interpolated rate where the median ratio crosses 1 (`inf` = never, within the tested range). `|dv|` is in minutes at a 3:00:00 marathon (`180 * delta log-time`).

**p\* (rate where identity error = sampling noise):** both = inf, break = inf, join = inf

| op | p (requested) | realised per-athlete rate | spearman(v, v0) | top-10% Jaccard | med abs(dv) min@3:00 | sens_ratio med | sens_ratio p75 | sens_ratio max |
|---|---|---|---|---|---|---|---|---|
| both | 0.001 | 0.0020 | 0.9999 | 1.000 | 0.004 | 0.034 | 0.043 | 0.122 |
| both | 0.002 | 0.0040 | 0.9998 | 1.000 | 0.008 | 0.048 | 0.058 | 0.121 |
| both | 0.005 | 0.0100 | 0.9995 | 1.000 | 0.020 | 0.076 | 0.089 | 0.146 |
| both | 0.010 | 0.0199 | 0.9991 | 1.000 | 0.035 | 0.109 | 0.125 | 0.185 |
| both | 0.020 | 0.0398 | 0.9982 | 0.926 | 0.064 | 0.156 | 0.175 | 0.315 |
| both | 0.050 | 0.0988 | 0.9957 | 0.926 | 0.134 | 0.237 | 0.269 | 0.412 |
| both | 0.100 | 0.1950 | 0.9900 | 0.926 | 0.234 | 0.345 | 0.385 | 0.543 |
| break | 0.001 | 0.0010 | 1.0000 | 1.000 | 0.002 | 0.020 | 0.027 | 0.060 |
| break | 0.002 | 0.0020 | 0.9999 | 1.000 | 0.003 | 0.029 | 0.037 | 0.085 |
| break | 0.005 | 0.0050 | 0.9998 | 1.000 | 0.009 | 0.047 | 0.055 | 0.091 |
| break | 0.010 | 0.0100 | 0.9997 | 1.000 | 0.016 | 0.065 | 0.074 | 0.173 |
| break | 0.020 | 0.0200 | 0.9994 | 1.000 | 0.030 | 0.095 | 0.108 | 0.176 |
| break | 0.050 | 0.0500 | 0.9985 | 0.926 | 0.061 | 0.149 | 0.169 | 0.254 |
| break | 0.100 | 0.1000 | 0.9969 | 0.926 | 0.102 | 0.211 | 0.240 | 0.326 |
| join | 0.001 | 0.0010 | 0.9999 | 1.000 | 0.003 | 0.027 | 0.034 | 0.108 |
| join | 0.002 | 0.0020 | 0.9999 | 1.000 | 0.005 | 0.037 | 0.046 | 0.115 |
| join | 0.005 | 0.0050 | 0.9997 | 1.000 | 0.013 | 0.061 | 0.073 | 0.128 |
| join | 0.010 | 0.0100 | 0.9994 | 1.000 | 0.025 | 0.084 | 0.099 | 0.174 |
| join | 0.020 | 0.0200 | 0.9989 | 1.000 | 0.049 | 0.120 | 0.140 | 0.262 |
| join | 0.050 | 0.0500 | 0.9969 | 0.926 | 0.111 | 0.189 | 0.223 | 0.348 |
| join | 0.100 | 0.1000 | 0.9928 | 0.926 | 0.200 | 0.269 | 0.311 | 0.491 |
