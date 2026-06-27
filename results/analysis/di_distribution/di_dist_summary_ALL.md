# d_i distribution summary (cohort ALL, n_i>=5)

*Only GAUGE-INVARIANT quantities are reportable: the level/mean of d_i is pinned by the EB prior, so shape (spread, skew) and slopes only.*

*d_i < 0 = improver; frac_improver and corr(d,u) are level-anchored and explicitly flagged NOT reportable.*

_Generated 2026-06-21._


## A. HETEROGENEITY MAGNITUDE  (gauge-invariant)

| Quantity | Value | Source |
|---|---|---|
| sqrt(omega_d2)=EB prior SD  [M] | 1.68 %/yr | `athlete_drift/omega_profile/ALL_M.../profile_nu8p00.parquet, is_free` |
| sqrt(omega_d2)=EB prior SD  [W] | 1.62 %/yr | `athlete_drift/omega_profile/ALL_W.../profile_nu8p00.parquet, is_free` |
| disatt SD (n>=5)  [M] | 1.29 [1.28, 1.30] %/yr | `di_descriptors.csv, disatt_sd(+CI)` |
| disatt SD (n>=5)  [W] | 1.23 [1.21, 1.25] %/yr | `di_descriptors.csv, disatt_sd(+CI)` |
| -> reading | a 1-SD athlete drifts ~1.3-1.7 %/yr off the shared curve (noise-corrected; gauge-free) | `derived` |

## B. SHAPE  (gauge-invariant)

| Quantity | Value | Source |
|---|---|---|
| skew (n>=5)  [M] | +0.103 [+0.068, +0.139]  (CI excludes 0) | `di_descriptors.csv, skew(+CI)` |
| skew (n>=5)  [W] | +0.120 [+0.069, +0.171]  (CI excludes 0) | `di_descriptors.csv, skew(+CI)` |
| -> reading | modestly right-skewed: a tail of rapid decliners (d>0); small in magnitude | `derived` |

## C. SEX COMPARISON  (only gauge-safe SD & skew)

| Quantity | Value | Source |
|---|---|---|
| disatt SD  M-W (n>=5) | +0.06 [+0.04, +0.08] %/yr | `di_descriptors.csv, group=M-W` |
| skew  M-W (n>=5) | -0.017 [-0.079, +0.047]  (CI includes 0) | `di_descriptors.csv, group=M-W` |
| -> reading | trajectory spread & shape essentially equal across sexes (no robust difference) | `derived` |

## D. ENTRY-AGE GRADIENT  (gauge-safe; slope is shift-invariant)

| Quantity | Value | Source |
|---|---|---|
| slope(d ~ entry age)  [M] | +0.11 [+0.09, +0.12] %/yr per 10yr   rho=+0.074 (var expl ~0.5%) | `di_entryage_slope.csv (ALL_M)` |
| slope(d ~ entry age)  [W] | +0.07 [+0.05, +0.09] %/yr per 10yr   rho=+0.052 (var expl ~0.3%) | `di_entryage_slope.csv (ALL_W)` |
| -> reading | gamma fan absorbs most of the entry-age effect; a faint residual remains (later debut -> decliner, Stones direction), an order of magnitude below the ~1.3 %/yr individual SD -> d_i is overwhelmingly idiosyncratic | `derived` |

## E. NOT REPORTABLE  (gauge/prior-anchored -- do NOT quote)

| Quantity | Value | Source |
|---|---|---|
| frac_improver | EXCLUDED: counts d<0, but the zero is prior-pinned (moves with the gauge); not comparable across slices | `caveat` |
| corr(d, u) | EXCLUDED: gauge-dependent via the u-side transform (u += -c*mean_i A_n); additionally era-confounded | `caveat` |
