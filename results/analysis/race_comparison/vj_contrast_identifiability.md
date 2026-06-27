# Cross-sex v_j contrast gauge bound (Dv = v_M - v_W on shared races)

Level and linear-in-date parts of Dv are PURE CONVENTION (G0/G1,
unpriced -- any value achievable); the table prices only the
quadratic (Gq) freedom against omega_d. `ratio` compares the R=1
plausible per-race quadratic gauge shift to the per-race bootstrap
SD of the contrast.

| Cohort | Mrc | n_shared | sd(Dv) | resid share | delta_kill | R_min | delta_plaus | shift_plaus med (min@3h) | boot sd med | ratio med | ratio p90 | Verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ALL | 2 | 337 | 0.0085 | 0.98 | -1.21e-04 | 0.05 | 2.65e-03 | 0.0230 (4.14) | 0.0037 | 5.03 | 14.87 | GAUGE > BOOT |
| ALL | 5 | 295 | 0.0089 | 0.96 | -5.86e-05 | 0.02 | 2.73e-03 | 0.0226 (4.08) | 0.0059 | 3.61 | 11.73 | GAUGE > BOOT |
| Po10 | 2 | 202 | 0.0092 | 0.78 | -1.07e-04 | 0.04 | 2.56e-03 | 0.0222 (3.99) | 0.0067 | 3.11 | 7.23 | GAUGE > BOOT |
| Po10 | 5 | 166 | 0.0100 | 0.75 | -1.01e-04 | 0.04 | 2.47e-03 | 0.0208 (3.74) | 0.0075 | 2.79 | 6.41 | GAUGE > BOOT |

_`sd(Dv)` = SD of the raw contrast over shared races (log units; x180 = minutes at 3:00). `resid share` = gauge-free fraction of the contrast variance (residual after projecting out {1, t, t^2}) -- the ceiling on any reportable cross-sex v_j pattern. `delta_kill` = quadratic-in-date coefficient that zeroes the fitted quadratic part; `R_min` = its drift cost as a fraction of the fitted drift spread, min over sexes (q02 convention: <0.5 fragile, >1.5 robust). `delta_plaus` = the R=1 quadratic re-pin via the cheapest sex; `shift_plaus med` = median per-race |delta_plaus * q(t)| it induces. `ratio` = that shift / per-race bootstrap SD of Dv: >1 means the priced gauge freedom alone exceeds sampling noise. The affine freedom is infinite regardless -- cross-sex v levels and date trends are never reportable._
