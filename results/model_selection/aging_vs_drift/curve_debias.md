# Aging-curve de-biasing: dedicated-runner restriction vs adding d_i

*On the full field the population aging curve is flattened by newcomers caught mid-improvement; d_i (full field) and the dedicated-runner cut are two independent corrections that restore the steeper common shape.*

*nu=8; mrc2 = everyone, mrc5 = dedicated runners (>=5 finishes).*

_Generated 2026-06-21._


## A. Peak (improvement) career-age, yr post-debut  (ALL_B_14-25_mrc2)

| Quantity | Value | Source |
|---|---|---|
| no-d, everyone (mrc2)  [the contaminated outlier] | 2.72 yr | `curve_debias.csv, peak_noD_mrc2` |
| +d, everyone (mrc2) | 2.48 yr | `curve_debias.csv, peak_withD_mrc2` |
| no-d, dedicated (mrc5) | 2.48 yr | `curve_debias.csv, peak_noD_mrc5` |
| +d, dedicated (mrc5) | 2.48 yr | `curve_debias.csv, peak_withD_mrc5` |
| -> +d on the full field reproduces the dedicated-runner peak | 2.48 yr  (vs dedicated 2.48 yr) | `derived` |

## B. Gap between everyone and dedicated-runner curves (RMS log-time)  (ALL_B_14-25_mrc2)

| Quantity | Value | Source |
|---|---|---|
| gap, no d_i | 0.0076 | `curve_debias.csv, gap_noD` |
| gap, with d_i | 0.0042 | `curve_debias.csv, gap_withD` |
| -> fraction of the no-d gap closed by d_i | 45% | `derived` |

## C. Cross-subset summary (gap no-d -> with-d, % closed)

| Quantity | Value | Source |
|---|---|---|
| ALL_B_14-25_mrc2 | 0.0076 -> 0.0042  (45% closed); +d peak 2.48 yr | `curve_debias.csv, gap_noD/gap_withD/peak_withD_mrc2` |
| ALL_M_14-25_mrc2 | 0.0085 -> 0.0052  (39% closed); +d peak 2.48 yr | `curve_debias.csv, gap_noD/gap_withD/peak_withD_mrc2` |
| ALL_W_14-25_mrc2 | 0.0057 -> 0.0020  (66% closed); +d peak 2.36 yr | `curve_debias.csv, gap_noD/gap_withD/peak_withD_mrc2` |
| Po10_B_14-25_mrc2 | 0.0062 -> 0.0023  (64% closed); +d peak 2.60 yr | `curve_debias.csv, gap_noD/gap_withD/peak_withD_mrc2` |
| Po10_M_14-25_mrc2 | 0.0055 -> 0.0026  (53% closed); +d peak 2.72 yr | `curve_debias.csv, gap_noD/gap_withD/peak_withD_mrc2` |
| Po10_W_14-25_mrc2 | 0.0052 -> 0.0023  (55% closed); +d peak 2.36 yr | `curve_debias.csv, gap_noD/gap_withD/peak_withD_mrc2` |
| WA_M_14-25_mrc2 | 0.0057 -> 0.0031  (45% closed); +d peak 2.47 yr | `curve_debias.csv, gap_noD/gap_withD/peak_withD_mrc2` |
| WA_W_14-25_mrc2 | 0.0065 -> 0.0025  (62% closed); +d peak 2.59 yr | `curve_debias.csv, gap_noD/gap_withD/peak_withD_mrc2` |
