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

slice: Po10_B_14-25_mrc2; min_n = 0; min_editions = 5; n_top = 20; model = full_nu8p00

### Po10_B_14-25_mrc2 -- fastest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | barcelona_marathon | ESP | 5 | 2021-25 | -0.0183 | [-0.0234, -0.0105] | -3.3 | 1.00 | [1, 10] | headline |
| 2 | frankfurt_marathon | GER | 9 | 2015-25 | -0.0172 | [-0.0211, -0.0108] | -3.1 | 1.00 | [1, 10] | headline |
| 3 | newport_marathon | GBR | 6 | 2018-25 | -0.0170 | [-0.0202, -0.0133] | -3.1 | 1.00 | [1, 8] | headline |
| 4 | dublin_marathon | IRL | 10 | 2014-25 | -0.0170 | [-0.0198, -0.0124] | -3.1 | 1.00 | [1, 9] | headline |
| 5 | seville_marathon | ESP | 7 | 2018-25 | -0.0166 | [-0.0202, -0.0114] | -3.0 | 1.00 | [1, 8] | headline |
| 6 | valencia_marathon | ESP | 8 | 2017-25 | -0.0158 | [-0.0171, -0.0127] | -2.8 | 1.00 | [2, 9] | headline |
| 7 | paris_marathon | FRA | 11 | 2014-25 | -0.0141 | [-0.0185, -0.0100] | -2.5 | 1.00 | [1, 12] | headline |
| 8 | manchester_marathon | GBR | 6 | 2019-25 | -0.0132 | [-0.0146, -0.0115] | -2.4 | 1.00 | [5, 10] | headline |
| 9 | bostonUK_marathon | GBR | 5 | 2021-25 | -0.0110 | [-0.0155, -0.0064] | -2.0 | 1.00 | [4, 16] | headline |
| 10 | amsterdam_marathon | NED | 9 | 2016-25 | -0.0106 | [-0.0164, -0.0082] | -1.9 | 1.00 | [4, 14] | headline |
| 11 | yorkshire_marathon | GBR | 6 | 2019-25 | -0.0103 | [-0.0125, -0.0062] | -1.9 | 1.00 | [8, 16] | headline |
| 12 | chester_marathon | GBR | 9 | 2016-25 | -0.0093 | [-0.0114, -0.0050] | -1.7 | 1.00 | [10, 16] | headline |
| 13 | berlin_marathon | GER | 11 | 2014-25 | -0.0091 | [-0.0102, -0.0066] | -1.6 | 1.00 | [9, 15] | headline |
| 14 | malaga_marathon | ESP | 7 | 2018-25 | -0.0073 | [-0.0137, -0.0009] | -1.3 | 0.99 | [6, 18] | headline |
| 15 | chicago_marathon | USA | 11 | 2014-25 | -0.0072 | [-0.0105, -0.0035] | -1.3 | 1.00 | [10, 17] | headline |
| 16 | rotterdam_marathon | NED | 6 | 2018-24 | -0.0038 | [-0.0127, +0.0023] | -0.7 | 0.97 | [9, 21] | headline |
| 17 | hamburg_marathon | GER | 6 | 2018-25 | -0.0031 | [-0.0136, +0.0074] | -0.5 | 0.85 | [6, 22] | headline |
| 18 | london_marathon | GBR | 11 | 2014-25 | -0.0010 | [-0.0035, -0.0006] | -0.2 | 1.00 | [16, 20] | headline |
| 19 | tokyo_marathon | JPN | 8 | 2015-25 | +0.0010 | [-0.0029, +0.0033] | +0.2 | 0.86 | [17, 22] | headline |
| 20 | belfast_marathon | GBR | 10 | 2014-25 | +0.0022 | [-0.0036, +0.0071] | +0.4 | 0.41 | [17, 23] | tie |
| 21 | edinburgh_marathon | GBR | 10 | 2014-25 | +0.0024 | [-0.0006, +0.0059] | +0.4 | 0.40 | [19, 22] | tie |
| 22 | copenhagen_marathon | DEN | 9 | 2015-25 | +0.0052 | [-0.0055, +0.0091] | +0.9 | 0.46 | [15, 23] | tie |

### Po10_B_14-25_mrc2 -- slowest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | boston_marathon | USA | 11 | 2014-25 | +0.0267 | [+0.0147, +0.0340] | +4.8 | 1.00 | [1, 4] | headline |
| 2 | nyc_marathon | USA | 11 | 2014-25 | +0.0205 | [+0.0176, +0.0238] | +3.7 | 1.00 | [1, 3] | headline |
| 3 | stockholm_marathon | SWE | 8 | 2014-24 | +0.0112 | [+0.0058, +0.0249] | +2.0 | 1.00 | [1, 5] | headline |
| 4 | vienna_marathon | AUT | 6 | 2018-25 | +0.0107 | [-0.0001, +0.0211] | +1.9 | 1.00 | [2, 9] | headline |
| 5 | milton_keynes_marathon | GBR | 11 | 2014-25 | +0.0081 | [+0.0044, +0.0152] | +1.5 | 1.00 | [3, 6] | headline |
| 6 | copenhagen_marathon | DEN | 9 | 2015-25 | +0.0052 | [-0.0055, +0.0091] | +0.9 | 1.00 | [5, 13] | headline |
| 7 | edinburgh_marathon | GBR | 10 | 2014-25 | +0.0024 | [-0.0006, +0.0059] | +0.4 | 1.00 | [6, 9] | headline |
| 8 | belfast_marathon | GBR | 10 | 2014-25 | +0.0022 | [-0.0036, +0.0071] | +0.4 | 1.00 | [5, 11] | headline |
| 9 | tokyo_marathon | JPN | 8 | 2015-25 | +0.0010 | [-0.0029, +0.0033] | +0.2 | 1.00 | [6, 11] | headline |
| 10 | london_marathon | GBR | 11 | 2014-25 | -0.0010 | [-0.0035, -0.0006] | -0.2 | 1.00 | [8, 12] | headline |
| 11 | hamburg_marathon | GER | 6 | 2018-25 | -0.0031 | [-0.0136, +0.0074] | -0.5 | 0.96 | [6, 22] | headline |
| 12 | rotterdam_marathon | NED | 6 | 2018-24 | -0.0038 | [-0.0127, +0.0023] | -0.7 | 1.00 | [7, 19] | headline |
| 13 | chicago_marathon | USA | 11 | 2014-25 | -0.0072 | [-0.0105, -0.0035] | -1.3 | 0.99 | [11, 18] | headline |
| 14 | malaga_marathon | ESP | 7 | 2018-25 | -0.0073 | [-0.0137, -0.0009] | -1.3 | 0.96 | [10, 22] | headline |
| 15 | berlin_marathon | GER | 11 | 2014-25 | -0.0091 | [-0.0102, -0.0066] | -1.6 | 1.00 | [13, 19] | headline |
| 16 | chester_marathon | GBR | 9 | 2016-25 | -0.0093 | [-0.0114, -0.0050] | -1.7 | 1.00 | [12, 18] | headline |
| 17 | yorkshire_marathon | GBR | 6 | 2019-25 | -0.0103 | [-0.0125, -0.0062] | -1.9 | 0.98 | [12, 20] | headline |
| 18 | amsterdam_marathon | NED | 9 | 2016-25 | -0.0106 | [-0.0164, -0.0082] | -1.9 | 0.78 | [14, 24] | headline |
| 19 | bostonUK_marathon | GBR | 5 | 2021-25 | -0.0110 | [-0.0155, -0.0064] | -2.0 | 0.89 | [12, 24] | headline |
| 20 | manchester_marathon | GBR | 6 | 2019-25 | -0.0132 | [-0.0146, -0.0115] | -2.4 | 0.63 | [18, 23] | headline |
| 21 | paris_marathon | FRA | 11 | 2014-25 | -0.0141 | [-0.0185, -0.0100] | -2.5 | 0.28 | [16, 27] | tie |
