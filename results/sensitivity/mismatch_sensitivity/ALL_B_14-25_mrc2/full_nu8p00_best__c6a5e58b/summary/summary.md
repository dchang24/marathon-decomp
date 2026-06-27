# Identity-mismatch sensitivity -- ALL_B_14-25_mrc2/full_nu8p00_best__c6a5e58b

Inject ADDITIONAL linkage error on top of the production resolution -- `op` = `break` (split one runner's records into two, mimicking a failed join), `join` (merge two different runners, mimicking a wrong match), or `both` -- at requested rate `p`, refit under Monte-Carlo replication, and compare the perturbed race factor `v_j` to the unperturbed baseline on the common race set (rank 1 = hardest).

**sens_ratio** = mismatch per-race SD / bootstrap per-race SD: the mismatch-induced movement relative to the athlete-sampling noise the analysis already carries. `sens_ratio < 1` means identity error perturbs `v_j` by LESS than ordinary sampling noise. **p\*** = interpolated rate where the median ratio crosses 1 (`inf` = never, within the tested range). `|dv|` is in minutes at a 3:00:00 marathon (`180 * delta log-time`).

**p\* (rate where identity error = sampling noise):** both = inf, break = inf, join = inf

| op | p (requested) | realised per-athlete rate | spearman(v, v0) | top-10% Jaccard | med abs(dv) min@3:00 | sens_ratio med | sens_ratio p75 | sens_ratio max |
|---|---|---|---|---|---|---|---|---|
| both | 0.001 | 0.0020 | 0.9998 | 1.000 | 0.009 | 0.072 | 0.116 | 1.453 |
| both | 0.002 | 0.0040 | 0.9995 | 1.000 | 0.016 | 0.092 | 0.133 | 1.398 |
| both | 0.005 | 0.0100 | 0.9989 | 0.944 | 0.027 | 0.146 | 0.206 | 2.302 |
| both | 0.010 | 0.0200 | 0.9970 | 0.944 | 0.042 | 0.209 | 0.292 | 3.289 |
| both | 0.020 | 0.0398 | 0.9970 | 0.944 | 0.058 | 0.271 | 0.346 | 3.374 |
| both | 0.050 | 0.0988 | 0.9933 | 0.944 | 0.105 | 0.399 | 0.497 | 2.847 |
| both | 0.100 | 0.1950 | 0.9864 | 0.842 | 0.174 | 0.521 | 0.639 | 2.155 |
| break | 0.001 | 0.0010 | 0.9999 | 1.000 | 0.004 | 0.026 | 0.030 | 0.128 |
| break | 0.002 | 0.0020 | 0.9999 | 1.000 | 0.005 | 0.034 | 0.039 | 0.096 |
| break | 0.005 | 0.0050 | 0.9999 | 1.000 | 0.010 | 0.057 | 0.063 | 0.189 |
| break | 0.010 | 0.0100 | 0.9998 | 1.000 | 0.015 | 0.079 | 0.090 | 0.141 |
| break | 0.020 | 0.0200 | 0.9997 | 1.000 | 0.022 | 0.111 | 0.123 | 0.180 |
| break | 0.050 | 0.0500 | 0.9994 | 0.944 | 0.039 | 0.177 | 0.200 | 0.309 |
| break | 0.100 | 0.1000 | 0.9988 | 0.944 | 0.062 | 0.254 | 0.281 | 0.398 |
| join | 0.001 | 0.0010 | 0.9998 | 1.000 | 0.009 | 0.060 | 0.091 | 1.062 |
| join | 0.002 | 0.0020 | 0.9996 | 1.000 | 0.013 | 0.083 | 0.125 | 1.366 |
| join | 0.005 | 0.0050 | 0.9987 | 0.944 | 0.025 | 0.133 | 0.196 | 2.135 |
| join | 0.010 | 0.0100 | 0.9968 | 0.944 | 0.038 | 0.188 | 0.274 | 3.262 |
| join | 0.020 | 0.0200 | 0.9972 | 0.944 | 0.051 | 0.242 | 0.318 | 3.293 |
| join | 0.050 | 0.0500 | 0.9935 | 0.944 | 0.093 | 0.353 | 0.441 | 2.792 |
| join | 0.100 | 0.1000 | 0.9870 | 0.842 | 0.159 | 0.461 | 0.591 | 2.250 |
