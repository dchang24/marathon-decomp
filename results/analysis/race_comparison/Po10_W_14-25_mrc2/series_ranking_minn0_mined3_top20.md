# Top-N fastest / slowest race series by median edition v_j

Estimand ('PB-chaser index'): per-series MEDIAN of edition v_j under
the beta=0 (bundling) gauge -- era-relative by construction, and
INCLUDING the series' typical weather (a net-of-weather course index
is deferred to the covariate analysis). The median discards one freak
edition per series with no hand-picked exclusions (hence the
min-editions floor of 3). Selection is by bootstrap rank stability
P(top-N), medians recomputed inside each replicate; no total order is
claimed -- ties are the P>=0.25 set. min@3:00 = 180*(exp(median_v)-1) = minutes
vs the average race for a 3:00:00 runner.

slice: Po10_W_14-25_mrc2; min_n = 0; min_editions = 3; n_top = 20; model = full_nu8p00

### Po10_W_14-25_mrc2 -- fastest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | amsterdam_marathon | NED | 8 | 2016-25 | -0.0181 | [-0.0220, -0.0102] | -3.2 | 1.00 | [1, 8] | headline |
| 2 | dublin_marathon | IRL | 10 | 2014-25 | -0.0177 | [-0.0209, -0.0109] | -3.2 | 1.00 | [1, 8] | headline |
| 3 | barcelona_marathon | ESP | 4 | 2021-25 | -0.0166 | [-0.0271, -0.0084] | -3.0 | 1.00 | [1, 11] | headline |
| 4 | newport_marathon | GBR | 6 | 2018-25 | -0.0155 | [-0.0209, -0.0076] | -2.8 | 1.00 | [1, 12] | headline |
| 5 | seville_marathon | ESP | 6 | 2019-25 | -0.0138 | [-0.0197, -0.0055] | -2.5 | 1.00 | [1, 13] | headline |
| 6 | frankfurt_marathon | GER | 7 | 2015-25 | -0.0132 | [-0.0230, -0.0102] | -2.4 | 1.00 | [1, 9] | headline |
| 7 | manchester_marathon | GBR | 6 | 2019-25 | -0.0117 | [-0.0131, -0.0079] | -2.1 | 1.00 | [5, 12] | headline |
| 8 | valencia_marathon | ESP | 8 | 2017-25 | -0.0110 | [-0.0152, -0.0057] | -2.0 | 1.00 | [3, 13] | headline |
| 9 | chester_marathon | GBR | 9 | 2016-25 | -0.0104 | [-0.0133, -0.0045] | -1.9 | 1.00 | [5, 14] | headline |
| 10 | paris_marathon | FRA | 10 | 2014-25 | -0.0099 | [-0.0170, -0.0018] | -1.8 | 1.00 | [3, 15] | headline |
| 11 | abingdon_marathon | GBR | 4 | 2022-25 | -0.0081 | [-0.0116, -0.0011] | -1.4 | 1.00 | [6, 16] | headline |
| 12 | yorkshire_marathon | GBR | 6 | 2019-25 | -0.0080 | [-0.0124, -0.0029] | -1.4 | 1.00 | [7, 15] | headline |
| 13 | malaga_marathon | ESP | 4 | 2018-25 | -0.0071 | [-0.0158, +0.0036] | -1.3 | 0.99 | [3, 19] | headline |
| 14 | bostonUK_marathon | GBR | 5 | 2021-25 | -0.0057 | [-0.0140, -0.0002] | -1.0 | 1.00 | [5, 18] | headline |
| 15 | berlin_marathon | GER | 11 | 2014-25 | -0.0039 | [-0.0066, +0.0004] | -0.7 | 1.00 | [12, 18] | headline |
| 16 | chicago_marathon | USA | 10 | 2014-25 | -0.0026 | [-0.0052, +0.0026] | -0.5 | 1.00 | [13, 19] | headline |
| 17 | london_marathon | GBR | 11 | 2014-25 | +0.0005 | [-0.0018, +0.0022] | +0.1 | 1.00 | [16, 19] | headline |
| 18 | tokyo_marathon | JPN | 7 | 2016-25 | +0.0009 | [-0.0034, +0.0068] | +0.2 | 0.95 | [14, 21] | headline |
| 19 | edinburgh_marathon | GBR | 10 | 2014-25 | +0.0031 | [-0.0004, +0.0074] | +0.6 | 0.93 | [16, 21] | headline |
| 20 | belfast_marathon | GBR | 10 | 2014-25 | +0.0031 | [-0.0018, +0.0114] | +0.6 | 0.82 | [16, 22] | headline |

### Po10_W_14-25_mrc2 -- slowest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | boston_marathon | USA | 10 | 2014-25 | +0.0226 | [+0.0169, +0.0282] | +4.1 | 1.00 | [1, 3] | headline |
| 2 | nyc_marathon | USA | 11 | 2014-25 | +0.0216 | [+0.0181, +0.0269] | +3.9 | 1.00 | [1, 3] | headline |
| 3 | copenhagen_marathon | DEN | 6 | 2018-25 | +0.0121 | [-0.0008, +0.0239] | +2.2 | 1.00 | [1, 9] | headline |
| 4 | brighton_marathon | GBR | 4 | 2022-25 | +0.0109 | [+0.0075, +0.0158] | +2.0 | 1.00 | [3, 5] | headline |
| 5 | milton_keynes_marathon | GBR | 11 | 2014-25 | +0.0063 | [+0.0040, +0.0176] | +1.1 | 1.00 | [3, 6] | headline |
| 6 | belfast_marathon | GBR | 10 | 2014-25 | +0.0031 | [-0.0018, +0.0114] | +0.6 | 1.00 | [4, 10] | headline |
| 7 | edinburgh_marathon | GBR | 10 | 2014-25 | +0.0031 | [-0.0004, +0.0074] | +0.6 | 1.00 | [5, 10] | headline |
| 8 | tokyo_marathon | JPN | 7 | 2016-25 | +0.0009 | [-0.0034, +0.0068] | +0.2 | 1.00 | [5, 12] | headline |
| 9 | london_marathon | GBR | 11 | 2014-25 | +0.0005 | [-0.0018, +0.0022] | +0.1 | 1.00 | [7, 10] | headline |
| 10 | chicago_marathon | USA | 10 | 2014-25 | -0.0026 | [-0.0052, +0.0026] | -0.5 | 1.00 | [7, 13] | headline |
| 11 | berlin_marathon | GER | 11 | 2014-25 | -0.0039 | [-0.0066, +0.0004] | -0.7 | 1.00 | [8, 14] | headline |
| 12 | bostonUK_marathon | GBR | 5 | 2021-25 | -0.0057 | [-0.0140, -0.0002] | -1.0 | 0.95 | [8, 21] | headline |
| 13 | malaga_marathon | ESP | 4 | 2018-25 | -0.0071 | [-0.0158, +0.0036] | -1.3 | 0.94 | [7, 23] | headline |
| 14 | yorkshire_marathon | GBR | 6 | 2019-25 | -0.0080 | [-0.0124, -0.0029] | -1.4 | 0.99 | [11, 19] | headline |
| 15 | abingdon_marathon | GBR | 4 | 2022-25 | -0.0081 | [-0.0116, -0.0011] | -1.4 | 0.98 | [10, 20] | headline |
| 16 | paris_marathon | FRA | 10 | 2014-25 | -0.0099 | [-0.0170, -0.0018] | -1.8 | 0.86 | [11, 23] | headline |
| 17 | chester_marathon | GBR | 9 | 2016-25 | -0.0104 | [-0.0133, -0.0045] | -1.9 | 0.94 | [12, 21] | headline |
| 18 | valencia_marathon | ESP | 8 | 2017-25 | -0.0110 | [-0.0152, -0.0057] | -2.0 | 0.81 | [13, 23] | headline |
| 19 | manchester_marathon | GBR | 6 | 2019-25 | -0.0117 | [-0.0131, -0.0079] | -2.1 | 0.93 | [14, 21] | headline |
| 20 | frankfurt_marathon | GER | 7 | 2015-25 | -0.0132 | [-0.0230, -0.0102] | -2.4 | 0.26 | [17, 25] | tie |
| 21 | seville_marathon | ESP | 6 | 2019-25 | -0.0138 | [-0.0197, -0.0055] | -2.5 | 0.47 | [13, 25] | tie |
| 22 | newport_marathon | GBR | 6 | 2018-25 | -0.0155 | [-0.0209, -0.0076] | -2.8 | 0.34 | [14, 25] | tie |
| 23 | barcelona_marathon | ESP | 4 | 2021-25 | -0.0166 | [-0.0271, -0.0084] | -3.0 | 0.25 | [15, 25] | tie |
