# Drift-prior variance omega_d2: identification and stability

*Per-athlete career drift d_i ~ N(0, omega_d2); omega_d2 is LEARNED (empirical Bayes), not chosen. nu=8, production operating point.*

*The marginal (not the plain data-fit) likelihood locates omega*; the EM update that finds it is a contraction to a single fixed point.*

_Generated 2026-06-21._


## A. EB-learned omega* across data subsets (stability)

| Quantity | Value | Source |
|---|---|---|
| ALL_B_14-25_mrc2 | omega* = 2.765e-04   (n_elig 171,990) | `omega_profile/ALL_B_14-25_mrc2/profile_nu8p00.parquet, is_free row, omega_d2` |
| ALL_M_14-25_mrc2 | omega* = 2.814e-04   (n_elig 120,919) | `omega_profile/ALL_M_14-25_mrc2/profile_nu8p00.parquet, is_free row, omega_d2` |
| ALL_W_14-25_mrc2 | omega* = 2.627e-04   (n_elig 51,040) | `omega_profile/ALL_W_14-25_mrc2/profile_nu8p00.parquet, is_free row, omega_d2` |
| Po10_B_14-25_mrc2 | omega* = 2.844e-04   (n_elig 27,396) | `omega_profile/Po10_B_14-25_mrc2/profile_nu8p00.parquet, is_free row, omega_d2` |
| Po10_M_14-25_mrc2 | omega* = 2.867e-04   (n_elig 18,140) | `omega_profile/Po10_M_14-25_mrc2/profile_nu8p00.parquet, is_free row, omega_d2` |
| Po10_W_14-25_mrc2 | omega* = 2.775e-04   (n_elig 9,123) | `omega_profile/Po10_W_14-25_mrc2/profile_nu8p00.parquet, is_free row, omega_d2` |
| range across subsets | 2.63e-04 - 2.87e-04  (spread 9%) | `derived` |

## B. Marginal likelihood locates omega*  (headline slice ALL_B_14-25_mrc2)

| Quantity | Value | Source |
|---|---|---|
| argmax_omega(logML) over the frozen sweep | omega_mult = 1  (i.e. at omega*) | `omega_profile/ALL_B_14-25_mrc2/profile_nu8p00.parquet, logML, argmax over frozen rows` |
| logML drop a decade BELOW omega* (x1/10) | 29,994 nats | `omega_profile/ALL_B_14-25_mrc2/profile_nu8p00.parquet, logML(x1) - logML(x0.1)` |
| logML drop a decade ABOVE omega* (x10) | 80,162 nats | `omega_profile/ALL_B_14-25_mrc2/profile_nu8p00.parquet, logML(x1) - logML(x10)` |
| plain data-fit log-lik monotone in omega? | yes  (so it never singles out omega*) | `omega_profile/ALL_B_14-25_mrc2/profile_nu8p00.parquet, data_loglik vs omega_mult` |

## C. Estimand insensitivity within omega* x [1/3, 3]  (ALL_B_14-25_mrc2)

| Quantity | Value | Source |
|---|---|---|
| max \|dv\| race factor over the band | 0.0073 log-time  (~0.7%) | `omega_profile/ALL_B_14-25_mrc2/profile_nu8p00.parquet, max_abs_dv, omega_mult in [1/3,3]` |
| max aging-curve deviation over the band | 0.0066 log-time  (~0.7%) | `omega_profile/ALL_B_14-25_mrc2/profile_nu8p00.parquet, aging_maxdev, omega_mult in [1/3,3]` |
| min race-factor correlation to the free fit over the band | 0.9951 | `omega_profile/ALL_B_14-25_mrc2/profile_nu8p00.parquet, corr_v_to_free, omega_mult in [1/3,3]` |

## D. omega* init-invariance (EM is a contraction)  (ALL_B_14-25_mrc2)

| Quantity | Value | Source |
|---|---|---|
| initial-guess range swept | 1e-07 - 1e-02  (6 inits) | `omega_init/ALL_B_14-25_mrc2/init_summary_nu8p00.parquet, omega_d2_init` |
| converged omega* spread (rel. to reference) | max \|omega/omega_ref - 1\| = 7.9e-04  over 6 converged inits | `omega_init/ALL_B_14-25_mrc2/init_summary_nu8p00.parquet, omega_final_rel_to_ref (a signed deviation, ref=0)` |
| iteration count range | 128 - 641 iters | `omega_init/ALL_B_14-25_mrc2/init_summary_nu8p00.parquet, n_iter (converged inits)` |
| inits that failed to converge | 0  (a too-small init can stall on the near-zero EM manifold; the failure is loud, never a silent wrong basin) | `omega_init/ALL_B_14-25_mrc2/init_summary_nu8p00.parquet, converged flag` |
