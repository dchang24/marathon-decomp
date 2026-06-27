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

slice: Po10_W_14-25_mrc2; min_n = 0; min_editions = 5; n_top = 20; model = full_nu8p00

### Po10_W_14-25_mrc2 -- fastest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | amsterdam_marathon | NED | 8 | 2016-25 | -0.0189 | [-0.0223, -0.0111] | -3.4 | 1.00 | [1, 5] | headline |
| 2 | newport_marathon | GBR | 6 | 2018-25 | -0.0157 | [-0.0206, -0.0092] | -2.8 | 1.00 | [1, 8] | headline |
| 3 | dublin_marathon | IRL | 10 | 2014-25 | -0.0154 | [-0.0211, -0.0106] | -2.8 | 1.00 | [1, 7] | headline |
| 4 | frankfurt_marathon | GER | 7 | 2015-25 | -0.0146 | [-0.0242, -0.0108] | -2.6 | 1.00 | [1, 7] | headline |
| 5 | seville_marathon | ESP | 6 | 2019-25 | -0.0125 | [-0.0191, -0.0051] | -2.2 | 1.00 | [1, 11] | headline |
| 6 | paris_marathon | FRA | 10 | 2014-25 | -0.0113 | [-0.0177, -0.0025] | -2.0 | 1.00 | [2, 13] | headline |
| 7 | valencia_marathon | ESP | 8 | 2017-25 | -0.0106 | [-0.0153, -0.0053] | -1.9 | 1.00 | [3, 11] | headline |
| 8 | manchester_marathon | GBR | 6 | 2019-25 | -0.0104 | [-0.0122, -0.0070] | -1.9 | 1.00 | [5, 10] | headline |
| 9 | chester_marathon | GBR | 9 | 2016-25 | -0.0100 | [-0.0137, -0.0045] | -1.8 | 1.00 | [4, 11] | headline |
| 10 | yorkshire_marathon | GBR | 6 | 2019-25 | -0.0062 | [-0.0122, -0.0023] | -1.1 | 1.00 | [6, 13] | headline |
| 11 | bostonUK_marathon | GBR | 5 | 2021-25 | -0.0048 | [-0.0139, +0.0010] | -0.9 | 1.00 | [4, 15] | headline |
| 12 | berlin_marathon | GER | 11 | 2014-25 | -0.0036 | [-0.0065, +0.0004] | -0.7 | 1.00 | [11, 14] | headline |
| 13 | tokyo_marathon | JPN | 7 | 2016-25 | -0.0009 | [-0.0052, +0.0052] | -0.2 | 1.00 | [11, 18] | headline |
| 14 | london_marathon | GBR | 11 | 2014-25 | -0.0009 | [-0.0022, +0.0019] | -0.2 | 1.00 | [12, 16] | headline |
| 15 | chicago_marathon | USA | 10 | 2014-25 | -0.0009 | [-0.0058, +0.0031] | -0.2 | 1.00 | [10, 17] | headline |
| 16 | belfast_marathon | GBR | 10 | 2014-25 | +0.0026 | [-0.0019, +0.0110] | +0.5 | 1.00 | [13, 19] | headline |
| 17 | edinburgh_marathon | GBR | 10 | 2014-25 | +0.0044 | [-0.0001, +0.0076] | +0.8 | 1.00 | [14, 18] | headline |
| 18 | milton_keynes_marathon | GBR | 11 | 2014-25 | +0.0065 | [+0.0015, +0.0170] | +1.2 | 0.99 | [15, 19] | headline |
| 19 | copenhagen_marathon | DEN | 6 | 2018-25 | +0.0119 | [-0.0006, +0.0239] | +2.1 | 0.95 | [14, 21] | headline |
| 20 | nyc_marathon | USA | 11 | 2014-25 | +0.0207 | [+0.0167, +0.0260] | +3.7 | 0.53 | [19, 21] | headline |
| 21 | boston_marathon | USA | 10 | 2014-25 | +0.0221 | [+0.0157, +0.0270] | +4.0 | 0.53 | [19, 21] | headline |

### Po10_W_14-25_mrc2 -- slowest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | boston_marathon | USA | 10 | 2014-25 | +0.0221 | [+0.0157, +0.0270] | +4.0 | 1.00 | [1, 3] | headline |
| 2 | nyc_marathon | USA | 11 | 2014-25 | +0.0207 | [+0.0167, +0.0260] | +3.7 | 1.00 | [1, 3] | headline |
| 3 | copenhagen_marathon | DEN | 6 | 2018-25 | +0.0119 | [-0.0006, +0.0239] | +2.1 | 1.00 | [1, 8] | headline |
| 4 | milton_keynes_marathon | GBR | 11 | 2014-25 | +0.0065 | [+0.0015, +0.0170] | +1.2 | 1.00 | [3, 7] | headline |
| 5 | edinburgh_marathon | GBR | 10 | 2014-25 | +0.0044 | [-0.0001, +0.0076] | +0.8 | 1.00 | [4, 8] | headline |
| 6 | belfast_marathon | GBR | 10 | 2014-25 | +0.0026 | [-0.0019, +0.0110] | +0.5 | 1.00 | [3, 9] | headline |
| 7 | chicago_marathon | USA | 10 | 2014-25 | -0.0009 | [-0.0058, +0.0031] | -0.2 | 1.00 | [5, 12] | headline |
| 8 | london_marathon | GBR | 11 | 2014-25 | -0.0009 | [-0.0022, +0.0019] | -0.2 | 1.00 | [6, 10] | headline |
| 9 | tokyo_marathon | JPN | 7 | 2016-25 | -0.0009 | [-0.0052, +0.0052] | -0.2 | 1.00 | [4, 11] | headline |
| 10 | berlin_marathon | GER | 11 | 2014-25 | -0.0036 | [-0.0065, +0.0004] | -0.7 | 1.00 | [8, 11] | headline |
| 11 | bostonUK_marathon | GBR | 5 | 2021-25 | -0.0048 | [-0.0139, +0.0010] | -0.9 | 1.00 | [7, 18] | headline |
| 12 | yorkshire_marathon | GBR | 6 | 2019-25 | -0.0062 | [-0.0122, -0.0023] | -1.1 | 1.00 | [9, 16] | headline |
| 13 | chester_marathon | GBR | 9 | 2016-25 | -0.0100 | [-0.0137, -0.0045] | -1.8 | 1.00 | [11, 18] | headline |
| 14 | manchester_marathon | GBR | 6 | 2019-25 | -0.0104 | [-0.0122, -0.0070] | -1.9 | 1.00 | [12, 17] | headline |
| 15 | valencia_marathon | ESP | 8 | 2017-25 | -0.0106 | [-0.0153, -0.0053] | -1.9 | 1.00 | [11, 19] | headline |
| 16 | paris_marathon | FRA | 10 | 2014-25 | -0.0113 | [-0.0177, -0.0025] | -2.0 | 0.98 | [9, 20] | headline |
| 17 | seville_marathon | ESP | 6 | 2019-25 | -0.0125 | [-0.0191, -0.0051] | -2.2 | 0.96 | [11, 21] | headline |
| 18 | frankfurt_marathon | GER | 7 | 2015-25 | -0.0146 | [-0.0242, -0.0108] | -2.6 | 0.70 | [15, 21] | headline |
| 19 | dublin_marathon | IRL | 10 | 2014-25 | -0.0154 | [-0.0211, -0.0106] | -2.8 | 0.84 | [15, 21] | headline |
| 20 | newport_marathon | GBR | 6 | 2018-25 | -0.0157 | [-0.0206, -0.0092] | -2.8 | 0.84 | [14, 21] | headline |
| 21 | amsterdam_marathon | NED | 8 | 2016-25 | -0.0189 | [-0.0223, -0.0111] | -3.4 | 0.68 | [17, 21] | headline |
