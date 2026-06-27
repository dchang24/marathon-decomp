# Top-N fastest / slowest race series by median edition v_j

Estimand ('PB-chaser index'): per-series MEDIAN of edition v_j under
the beta=0 (bundling) gauge -- era-relative by construction, and
INCLUDING the series' typical weather (a net-of-weather course index
is deferred to the covariate analysis). The median discards one freak
edition per series with no hand-picked exclusions (hence the
min-editions floor of 3). Selection is by bootstrap rank stability
P(top-N), medians recomputed inside each replicate; no total order is
claimed -- ties are the P>=0.25 set. min@3:00 = 180*median_v = minutes
vs the average race for a 3:00:00 runner.

slice: Po10_M_14-25_mrc2; min_n = 0; min_editions = 5; n_top = 20; model = full_nu8p00

### Po10_M_14-25_mrc2 -- fastest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | seville_marathon | ESP | 7 | 2018-25 | -0.0205 | [-0.0235, -0.0136] | -3.7 | 1.00 | [1, 7] | headline |
| 2 | newport_marathon | GBR | 6 | 2018-25 | -0.0180 | [-0.0208, -0.0132] | -3.2 | 1.00 | [1, 9] | headline |
| 3 | valencia_marathon | ESP | 8 | 2017-25 | -0.0178 | [-0.0202, -0.0128] | -3.2 | 1.00 | [1, 8] | headline |
| 4 | amsterdam_marathon | NED | 9 | 2016-25 | -0.0160 | [-0.0188, -0.0067] | -2.9 | 1.00 | [1, 14] | headline |
| 5 | dublin_marathon | IRL | 10 | 2014-25 | -0.0158 | [-0.0192, -0.0119] | -2.8 | 1.00 | [2, 10] | headline |
| 6 | manchester_marathon | GBR | 6 | 2019-25 | -0.0149 | [-0.0166, -0.0128] | -2.7 | 1.00 | [3, 9] | headline |
| 7 | frankfurt_marathon | GER | 9 | 2015-25 | -0.0149 | [-0.0198, -0.0111] | -2.7 | 1.00 | [1, 10] | headline |
| 8 | yorkshire_marathon | GBR | 6 | 2019-25 | -0.0137 | [-0.0162, -0.0089] | -2.5 | 1.00 | [4, 13] | headline |
| 9 | bostonUK_marathon | GBR | 5 | 2021-25 | -0.0119 | [-0.0196, -0.0081] | -2.1 | 1.00 | [1, 13] | headline |
| 10 | hamburg_marathon | GER | 5 | 2018-25 | -0.0116 | [-0.0216, +0.0084] | -2.1 | 0.95 | [1, 22] | headline |
| 11 | berlin_marathon | GER | 11 | 2014-25 | -0.0114 | [-0.0130, -0.0089] | -2.0 | 1.00 | [8, 13] | headline |
| 12 | chester_marathon | GBR | 9 | 2016-25 | -0.0104 | [-0.0116, -0.0040] | -1.9 | 1.00 | [9, 16] | headline |
| 13 | paris_marathon | FRA | 11 | 2014-25 | -0.0098 | [-0.0166, -0.0056] | -1.8 | 1.00 | [3, 15] | headline |
| 14 | chicago_marathon | USA | 11 | 2014-25 | -0.0073 | [-0.0123, -0.0041] | -1.3 | 1.00 | [9, 15] | headline |
| 15 | malaga_marathon | ESP | 5 | 2018-25 | -0.0068 | [-0.0136, +0.0013] | -1.2 | 1.00 | [8, 20] | headline |
| 16 | london_marathon | GBR | 11 | 2014-25 | -0.0027 | [-0.0046, -0.0012] | -0.5 | 1.00 | [14, 18] | headline |
| 17 | edinburgh_marathon | GBR | 10 | 2014-25 | +0.0000 | [-0.0038, +0.0030] | +0.0 | 0.94 | [15, 21] | headline |
| 18 | rotterdam_marathon | NED | 5 | 2018-24 | +0.0002 | [-0.0087, +0.0090] | +0.0 | 0.72 | [13, 22] | headline |
| 19 | copenhagen_marathon | DEN | 9 | 2015-25 | +0.0006 | [-0.0084, +0.0087] | +0.1 | 0.86 | [12, 22] | headline |
| 20 | vienna_marathon | AUT | 6 | 2018-25 | +0.0025 | [-0.0073, +0.0189] | +0.4 | 0.46 | [13, 24] | tie |
| 21 | tokyo_marathon | JPN | 8 | 2015-25 | +0.0026 | [-0.0027, +0.0059] | +0.5 | 0.80 | [16, 22] | headline |
| 22 | belfast_marathon | GBR | 10 | 2014-25 | +0.0063 | [-0.0023, +0.0103] | +1.1 | 0.26 | [16, 22] | tie |

### Po10_M_14-25_mrc2 -- slowest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | boston_marathon | USA | 11 | 2014-25 | +0.0309 | [+0.0131, +0.0339] | +5.6 | 1.00 | [1, 4] | headline |
| 2 | nyc_marathon | USA | 11 | 2014-25 | +0.0226 | [+0.0178, +0.0255] | +4.1 | 1.00 | [1, 3] | headline |
| 3 | milton_keynes_marathon | GBR | 11 | 2014-25 | +0.0155 | [+0.0087, +0.0189] | +2.8 | 1.00 | [2, 5] | headline |
| 4 | stockholm_marathon | SWE | 8 | 2014-24 | +0.0134 | [+0.0044, +0.0289] | +2.4 | 1.00 | [1, 6] | headline |
| 5 | belfast_marathon | GBR | 10 | 2014-25 | +0.0063 | [-0.0023, +0.0103] | +1.1 | 1.00 | [5, 11] | headline |
| 6 | tokyo_marathon | JPN | 8 | 2015-25 | +0.0026 | [-0.0027, +0.0059] | +0.5 | 1.00 | [5, 11] | headline |
| 7 | vienna_marathon | AUT | 6 | 2018-25 | +0.0025 | [-0.0073, +0.0189] | +0.4 | 1.00 | [3, 14] | headline |
| 8 | copenhagen_marathon | DEN | 9 | 2015-25 | +0.0006 | [-0.0084, +0.0087] | +0.1 | 1.00 | [5, 15] | headline |
| 9 | rotterdam_marathon | NED | 5 | 2018-24 | +0.0002 | [-0.0087, +0.0090] | +0.0 | 0.99 | [5, 14] | headline |
| 10 | edinburgh_marathon | GBR | 10 | 2014-25 | +0.0000 | [-0.0038, +0.0030] | +0.0 | 1.00 | [6, 12] | headline |
| 11 | london_marathon | GBR | 11 | 2014-25 | -0.0027 | [-0.0046, -0.0012] | -0.5 | 1.00 | [9, 13] | headline |
| 12 | malaga_marathon | ESP | 5 | 2018-25 | -0.0068 | [-0.0136, +0.0013] | -1.2 | 0.99 | [7, 19] | headline |
| 13 | chicago_marathon | USA | 11 | 2014-25 | -0.0073 | [-0.0123, -0.0041] | -1.3 | 1.00 | [12, 18] | headline |
| 14 | paris_marathon | FRA | 11 | 2014-25 | -0.0098 | [-0.0166, -0.0056] | -1.8 | 0.82 | [12, 24] | headline |
| 15 | chester_marathon | GBR | 9 | 2016-25 | -0.0104 | [-0.0116, -0.0040] | -1.9 | 1.00 | [11, 18] | headline |
| 16 | berlin_marathon | GER | 11 | 2014-25 | -0.0114 | [-0.0130, -0.0089] | -2.0 | 1.00 | [14, 19] | headline |
| 17 | hamburg_marathon | GER | 5 | 2018-25 | -0.0116 | [-0.0216, +0.0084] | -2.1 | 0.78 | [5, 26] | headline |
| 18 | bostonUK_marathon | GBR | 5 | 2021-25 | -0.0119 | [-0.0196, -0.0081] | -2.1 | 0.63 | [14, 26] | headline |
| 19 | yorkshire_marathon | GBR | 6 | 2019-25 | -0.0137 | [-0.0162, -0.0089] | -2.5 | 0.88 | [14, 23] | headline |
| 20 | frankfurt_marathon | GER | 9 | 2015-25 | -0.0149 | [-0.0198, -0.0111] | -2.7 | 0.38 | [17, 26] | tie |
| 21 | manchester_marathon | GBR | 6 | 2019-25 | -0.0149 | [-0.0166, -0.0128] | -2.7 | 0.36 | [18, 24] | tie |
| 22 | dublin_marathon | IRL | 10 | 2014-25 | -0.0158 | [-0.0192, -0.0119] | -2.8 | 0.34 | [17, 25] | tie |
| 23 | amsterdam_marathon | NED | 9 | 2016-25 | -0.0160 | [-0.0188, -0.0067] | -2.9 | 0.57 | [13, 26] | headline |
