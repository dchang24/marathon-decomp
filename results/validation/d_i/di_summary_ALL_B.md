# d_i inclusion / v_j de-biasing summary (ALL_B, nu=8)

*The per-athlete career drift d_i is justified by DE-BIASING v_j (it removes the field's career-stage composition leak), not by predictive fit -- so the headline is the v_j-bias test, not an information criterion.*

*Reads only the d_i validation + aging-vs-drift + athlete-drift outputs; never refits.*

_Generated 2026-06-21._


## A. V_J DE-BIASING  (headline: EB, everyone-field mrc2)

| Quantity | Value | Source |
|---|---|---|
| slug / n_races | ALL_B_14-25_mrc2 / 347 | `bias_corr.csv` |
| year-partialled corr(delta,v)  no-d | +0.2038 | `bias_summary.csv, eb, corr_nod` |
| +d | +0.0643 | `bias_summary.csv, eb, corr_d` |
| -> bias removed | 68%  (+0.204 -> +0.064) | `bias_summary.csv, eb, frac_removed` |
| partial slope  v ~ delta   no-d | 1.502 | `bias_corr.csv, eb/no-d, slope_partial` |
| +d | 0.461 | `bias_corr.csv, eb/+d, slope_partial  (1:1 absorption -> collapses)` |
| raw vs year-partialled corr  no-d | +0.075 -> +0.204 | `bias_corr.csv, eb/no-d, pearson_raw / pearson_partial` |
| (why partial>raw) | era trend in v + window trend in delta cancel in raw -> partialling exposes the coupling | `derived` |

## B. INDEPENDENT CONFIRMATION  (LOO + permutation, clean mrc5)

| Quantity | Value | Source |
|---|---|---|
| LOO year-partialled corr  no-d / +d | +0.2318 / +0.0726 | `bias_summary.csv, loo (mrc5), corr_nod/corr_d` |
| -> bias removed (LOO) | 69% | `bias_summary.csv, loo, frac_removed` |
| LOO partial slope  no-d / +d | 1.187 / 0.355 | `bias_corr.csv, loo (mrc5), slope_partial` |
| estimator agreement corr(EB,LOO)  mrc5 / mrc2 | 0.724 / 0.396 | `delta_agreement.csv, corr_eb_loo` |
| permutation null (mrc5)  no-d | z=4.81  p=0.0020 | `permutation/ALL_B.csv, no-d, z/p_one` |
| +d | z=1.73  p=0.0579  (loses signif.) | `permutation/ALL_B.csv, +d, z/p_one` |
| permutation pool size | 500 draws  (p floor = 1/(N+1)) | `permutation/ALL_B.csv, n_perm_total` |

## C. CROSS-CELL CONSISTENCY  (all Po10/ALL EB cells)

| Quantity | Value | Source |
|---|---|---|
| EB no-d corr range over cells | +0.186 .. +0.311  (n=12 cells) | `bias_summary.csv, all eb rows, corr_nod` |
| min / max cell | ALL_W_14-25_mrc2 / ALL_W_14-25_mrc5 | `derived` |
| EB frac_removed range | 0.61 .. 1.09 | `bias_summary.csv, all eb rows, frac_removed  (>1 = mrc2 LOO artifact n/a here)` |

## D. AGING-CURVE CORROBORATION  (cross-mrc de-biasing)

| Quantity | Value | Source |
|---|---|---|
| everyone-vs-dedicated curve gap  no-d | 0.0076 | `curve_debias.csv, gap_noD  (log-time RMS)` |
| +d | 0.0042 | `curve_debias.csv, gap_withD` |
| -> gap closed by d_i | 45% | `derived` |
| d_i curve shift  everyone(mrc2) / dedicated(mrc5) | 0.0074 / 0.0041  (more work where contaminated) | `curve_debias.csv, d_on_mrc2 / d_on_mrc5` |
| peak age (yr post-debut)  noD-mrc2 / +d-mrc2 / noD-mrc5 | 2.72 / 2.48 / 2.48  (+d-mrc2 == noD-mrc5) | `curve_debias.csv, peak_* (to grid resolution)` |
| curvature preserved (aging vs full, mrc2) | curve_corr 0.9948 | `block_shift.csv, curve_corr  (shape kept; only tilt removed)` |
| n eligible athletes with free d_i (mrc2) | 172,127 | `block_shift.csv, n_eligible` |

## E. PRIOR IDENTIFICATION  (omega_d2; numerical 'it works')

| Quantity | Value | Source |
|---|---|---|
| EB-learned omega_d2*  (mrc2) | 2.765e-04 | `omega_profile/ALL_B_14-25_mrc2/profile_nu8p00.parquet, is_free row` |
| marginal logML peak at | omega_mult = 1  (1.0 = omega*) | `argmax logML  (peaked, not flat -> identified)` |
| n eligible (free d_i) | 171,990 | `profile, n_elig` |
| insensitivity in omega* x [1/3,3] | aging dev <= 0.0066,  corr_v_to_free >= 0.9951 | `profile, aging_maxdev / corr_v_to_free  (<1% move)` |
| init-invariance (converged inits) | omega* rel-spread 9.3e-04,  aging dev <= 3.6e-05 | `omega_init/ALL_B_14-25_mrc2/init_summary_nu8p00.parquet, converged rows` |
| => d_i adds no tuning knob | prior learned (EB), identified, init-invariant | `derived` |
