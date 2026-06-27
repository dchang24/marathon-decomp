# Single-nu decision for the rank-1 baseline

*held-out CV/cell (1-SE rule) + v-stability plateau across data subsets; no new fitting (reads the e01_nu_cv CV sweep).*

_Generated 2026-06-21._


## Per-slice nu (1-SE rule on held-out CV/cell)

| Quantity | Value | Source |
|---|---|---|
| ALL_B_14-25_mrc2 | argmax 8; 1-SE band [8, 10]; CV* 0.9077 +- 0.0012 | `ALL_B_14-25_mrc2/cv_folds.csv (heldout_per_cell, K folds)` |
| ALL_M_14-25_mrc2 | argmax 8; 1-SE band [8, 10]; CV* 0.9002 +- 0.0036 | `ALL_M_14-25_mrc2/cv_folds.csv (heldout_per_cell, K folds)` |
| ALL_W_14-25_mrc2 | argmax 8; 1-SE band [6, 10]; CV* 0.9222 +- 0.0048 | `ALL_W_14-25_mrc2/cv_folds.csv (heldout_per_cell, K folds)` |
| Po10_B_14-25_mrc2 | argmax 6; 1-SE band [5, 8]; CV* 0.8873 +- 0.0040 | `Po10_B_14-25_mrc2/cv_folds.csv (heldout_per_cell, K folds)` |
| Po10_M_14-25_mrc2 | argmax 6; 1-SE band [5, 6]; CV* 0.8541 +- 0.0033 | `Po10_M_14-25_mrc2/cv_folds.csv (heldout_per_cell, K folds)` |
| Po10_W_14-25_mrc2 | argmax 6; 1-SE band [5, 8]; CV* 0.9453 +- 0.0088 | `Po10_W_14-25_mrc2/cv_folds.csv (heldout_per_cell, K folds)` |
| WA_M_14-25_mrc2 | argmax 5; 1-SE band [3, 10]; CV* 1.6155 +- 0.0239 | `WA_M_14-25_mrc2/cv_folds.csv (heldout_per_cell, K folds)` |
| WA_W_14-25_mrc2 | argmax 8; 1-SE band [4, 10]; CV* 1.4115 +- 0.0177 | `WA_W_14-25_mrc2/cv_folds.csv (heldout_per_cell, K folds)` |

## Recommended single shared nu

| Quantity | Value | Source |
|---|---|---|
| recommended nu* | 8 | `derived: grid value in the most 1-SE bands` |
| 1-SE interval intersection (all slices) | empty | `derived from per-slice 1-SE bands` |
| coverage at nu* | 7/8 slices' 1-SE bands contain nu* | `derived` |
| max deliverable v-cost at nu* (vs each slice's argmax-v) | 1-Pearson 2.11e-03, 1-Spearman 2.84e-03, max\|dv\| 6.37e-03 (log-time) | `derived from v_xnu.csv` |

## Coverage by grid nu (# slices whose 1-SE band contains it)

| Quantity | Value | Source |
|---|---|---|
| nu = 3 | 1/8 | `derived` |
| nu = 4 | 2/8 | `derived` |
| nu = 5 | 5/8 | `derived` |
| nu = 6 | 6/8 | `derived` |
| nu = 8 | 7/8 | `derived` |
| nu = 10 | 5/8 | `derived` |
| nu = 15 | 0/8 | `derived` |

## Notes

- The single shared nu rests on the v-plateau when the per-slice 1-SE intervals do not all overlap: CV can resolve nearby nu, but the deliverable race factor v barely moves (see max v-cost above).
- Folds are stratified within-athlete, so the across-fold SE slightly *over*states the spread (a principled heuristic, not a strict CI).
