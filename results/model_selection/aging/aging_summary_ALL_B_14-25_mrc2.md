# Aging progression form selection (ALL_B_14-25_mrc2)

*Production form = spline4-gvarying (natural cubic spline, 4 knots, varying gamma); nu=8; drift d_i OFF during selection.*

*Decision: held-out CV/cell primary -> curve plausibility guardrail -> in-sample BIC tie-break.*

_Generated 2026-06-21._


## A. AGING BLOCK vs NONE  (does an aging term help at all?)

| Quantity | Value | Source |
|---|---|---|
| held-out logdens/cell  (no aging) | 0.9077 | `form_compare.csv (slug, nu=8), form='[baseline no-aging]', cv_per_cell` |
| held-out logdens/cell  (spline4-gvarying) | 0.9263 | `form_compare.csv (slug, nu=8), form='spline4-gvarying', cv_per_cell` |
| -> CV gain from adding aging | +0.0186 nats/cell | `derived` |
| paired z (no-aging vs CV-best, same folds) | 23.9 | `form_compare.csv (slug, nu=8), form='[baseline no-aging]', paired_z` |
| in-sample BIC  (no aging / chosen) | 1,825,047 / 1,766,038 | `form_compare.csv (slug, nu=8), bic` |
| => term earns its keep? | in-sample BIC -59,009 AND held-out CV +0.0186  (both improve) | `derived` |

## B. WHICH FORM  (differences among flexible forms are tiny)

| Quantity | Value | Source |
|---|---|---|
| CV-best form on this slice | spline6-gvarying  (CV 0.9269) | `best_form_all.csv (slug, nu=8), best_cand / cv_per_cell` |
| chosen spline4-gvarying CV/cell | 0.9263 | `form_compare.csv (slug, nu=8), form='spline4-gvarying', cv_per_cell` |
| -> CV shortfall vs CV-best | 5.40e-04 nats/cell  (z 5.0) | `form_compare.csv (slug, nu=8), dCV_vs_best / paired_z` |
| AIC / BIC  spline4-gvarying | -3,083,081 / 1,766,038 | `form_compare.csv (slug, nu=8), form='spline4-gvarying', aic / bic` |
| effective d.o.f. / naive n_params | 402,949 / 402,950  (equal: no EB shrinkage, drift off) | `derived / form_compare.csv (slug, nu=8), k` |
| BIC  CV-best (spline6-gvarying) | 1,765,482 | `form_compare.csv (slug, nu=8), form='spline6-gvarying', bic` |
| -> BIC gap chosen vs CV-best | 555  (scale ~1.8e6) | `derived` |

## C. CURVE PLAUSIBILITY GUARDRAIL  (spline vs high-order polynomial)

| Quantity | Value | Source |
|---|---|---|
| tail non-monotone?  poly5-gvarying | 1  (sign changes 1) | `curve_metrics_all.csv (slug, nu=8, entry_age=mean), cand='poly5-gvarying', tail_nonmonotone / tail_sign_changes` |
| tail non-monotone?  spline4-gvarying | 0  (sign changes 0) | `curve_metrics_all.csv (slug, nu=8, entry_age=mean), cand='spline4-gvarying', tail_nonmonotone / tail_sign_changes` |
| peak (improvement) age  spline4-gvarying | 42.1 yr | `curve_metrics_all.csv (slug, nu=8, entry_age=mean), cand='spline4-gvarying', peak_age` |

## D. ENTRY-AGE INTERACTION  (varying gamma, after Stones 2019)

| Quantity | Value | Source |
|---|---|---|
| BIC  no entry-age term (spline4-goff) | 1,804,923 | `form_compare.csv (slug, nu=8), form='spline4-goff', bic` |
| BIC  varying gamma (spline4-gvarying) | 1,766,038 | `form_compare.csv (slug, nu=8), form='spline4-gvarying', bic` |
| -> BIC improvement from entry-age term | 38,885 | `derived` |
| CV/cell  scalar gamma (spline4-gscalar) | 0.9231  (z 14.5) | `form_compare.csv (slug, nu=8), form='spline4-gscalar', cv_per_cell / paired_z` |
| CV/cell  varying gamma (spline4-gvarying) | 0.9263 | `form_compare.csv (slug, nu=8), form='spline4-gvarying', cv_per_cell` |

## E. v ROBUST TO FORM  (the target output barely moves)

| Quantity | Value | Source |
|---|---|---|
| v Pearson   spline4-gvarying vs poly2-gscalar | 0.9960 | `v_vs_reference_all.csv (slug, nu=8, cand='spline4-gvarying', ref='poly2-gscalar'), pearson` |
| v Spearman  spline4-gvarying vs poly2-gscalar | 0.9944 | `v_vs_reference_all.csv (slug, nu=8, cand='spline4-gvarying', ref='poly2-gscalar'), spearman` |
| mean \|dv\| vs reference (log-time) | 0.0013 | `v_vs_reference_all.csv (slug, nu=8, cand='spline4-gvarying', ref='poly2-gscalar'), mean_abs_dv` |

## F. WHY spline4 (vs spline3 / spline5,6 / poly4,6)

| Quantity | Value | Source |
|---|---|---|
| spline3-gvarying | dof 402,947  BIC 1,770,164  CV 0.925078  capt  90.7% | `grid/metrics.csv + cv/form_selection.csv` |
| spline4-gvarying | dof 402,949  BIC 1,766,038  CV 0.926330  capt  97.2%  <- CHOSEN | `grid/metrics.csv + cv/form_selection.csv` |
| spline5-gvarying | dof 402,951  BIC 1,765,532  CV 0.926509  capt  98.1% | `grid/metrics.csv + cv/form_selection.csv` |
| spline6-gvarying | dof 402,953  BIC 1,765,482  CV 0.926870  capt 100.0%  <- CV-best | `grid/metrics.csv + cv/form_selection.csv` |
| poly4-gvarying | dof 402,951  BIC 1,766,461  CV   n/a  (not CV-scored) | `grid/metrics.csv + cv/form_selection.csv` |
| poly6-gvarying | dof 402,955  BIC 1,765,593  CV   n/a  (not CV-scored) | `grid/metrics.csv + cv/form_selection.csv` |
| spline4 vs spline3 (too simple) | CV +0.0013/cell (91%->97% capt), BIC -4,126 | `derived` |
| spline4 vs spline6 (CV-best) | CV 5.40e-04/cell (last 2.8%), BIC +555 | `derived` |
| spline4 vs poly4 / poly6 (in-sample BIC only) | dBIC -423 / +445  (<0 = spline4 better; poly route rejected on tail shape, see C) | `derived` |
| v agreement spline4-gv vs spline6-gv | Pearson 0.9997, Spearman 0.9995 | `ALL_B_14-25_mrc2/v_xform_agreement.csv` |
| CV-best form across ALL slices | ALL_B:spline6-gvarying, ALL_M:spline5-gvarying, ALL_W:spline6-gvarying | `best_form_all.csv (nu=8, ALL_*) -- unstable knot count, all within ~5e-4 of spline4` |
