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

slice: ALL_B_14-25_mrc2; min_n = 0; min_editions = 5; n_top = 20; model = full_nu8p00

### ALL_B_14-25_mrc2 -- fastest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | newport_marathon | GBR | 6 | 2018-25 | -0.0169 | [-0.0199, -0.0129] | -3.0 | 1.00 | [1, 5] | headline |
| 2 | frankfurt_marathon | GER | 9 | 2015-25 | -0.0153 | [-0.0169, -0.0137] | -2.8 | 1.00 | [1, 5] | headline |
| 3 | amsterdam_marathon | NED | 9 | 2016-25 | -0.0133 | [-0.0151, -0.0114] | -2.4 | 1.00 | [2, 7] | headline |
| 4 | eindhoven_marathon | NED | 8 | 2017-25 | -0.0128 | [-0.0159, -0.0104] | -2.3 | 1.00 | [1, 10] | headline |
| 5 | seville_marathon | ESP | 7 | 2018-25 | -0.0124 | [-0.0157, -0.0101] | -2.2 | 1.00 | [2, 10] | headline |
| 6 | dublin_marathon | IRL | 10 | 2014-25 | -0.0121 | [-0.0171, -0.0096] | -2.2 | 1.00 | [1, 11] | headline |
| 7 | zurich_marathon | SUI | 5 | 2019-25 | -0.0121 | [-0.0147, -0.0077] | -2.2 | 1.00 | [2, 14] | headline |
| 8 | hannover_marathon | GER | 6 | 2018-25 | -0.0120 | [-0.0140, -0.0095] | -2.2 | 1.00 | [4, 12] | headline |
| 9 | paris_marathon | FRA | 11 | 2014-25 | -0.0106 | [-0.0131, -0.0091] | -1.9 | 1.00 | [5, 13] | headline |
| 10 | manchester_marathon | GBR | 6 | 2019-25 | -0.0098 | [-0.0115, -0.0082] | -1.8 | 1.00 | [7, 14] | headline |
| 11 | valencia_marathon | ESP | 9 | 2017-25 | -0.0097 | [-0.0115, -0.0085] | -1.7 | 1.00 | [7, 14] | headline |
| 12 | berlin_marathon | GER | 11 | 2014-25 | -0.0084 | [-0.0090, -0.0067] | -1.5 | 1.00 | [12, 19] | headline |
| 13 | malaga_marathon | ESP | 7 | 2018-25 | -0.0083 | [-0.0112, -0.0025] | -1.5 | 0.76 | [8, 24] | headline |
| 14 | bostonUK_marathon | GBR | 5 | 2021-25 | -0.0077 | [-0.0132, -0.0051] | -1.4 | 0.94 | [5, 22] | headline |
| 15 | chicago_marathon | USA | 11 | 2014-25 | -0.0077 | [-0.0090, -0.0067] | -1.4 | 1.00 | [12, 19] | headline |
| 16 | barcelona_marathon | ESP | 5 | 2021-25 | -0.0076 | [-0.0106, -0.0047] | -1.4 | 0.87 | [9, 22] | headline |
| 17 | rotterdam_marathon | NED | 7 | 2018-25 | -0.0073 | [-0.0087, -0.0056] | -1.3 | 0.98 | [12, 20] | headline |
| 18 | stockholm_marathon | SWE | 11 | 2014-25 | -0.0070 | [-0.0095, -0.0049] | -1.3 | 0.92 | [11, 22] | headline |
| 19 | hamburg_marathon | GER | 7 | 2018-25 | -0.0064 | [-0.0070, -0.0042] | -1.1 | 0.54 | [18, 23] | headline |
| 20 | tokyo_marathon | JPN | 9 | 2015-25 | -0.0061 | [-0.0084, -0.0046] | -1.1 | 0.82 | [13, 22] | headline |
| 21 | yorkshire_marathon | GBR | 6 | 2019-25 | -0.0057 | [-0.0084, -0.0026] | -1.0 | 0.48 | [13, 24] | tie |
| 22 | prague_marathon | CZE | 10 | 2014-25 | -0.0053 | [-0.0098, +0.0000] | -0.9 | 0.43 | [10, 25] | tie |
| 23 | chester_marathon | GBR | 9 | 2016-25 | -0.0047 | [-0.0072, -0.0020] | -0.8 | 0.26 | [17, 24] | tie |

### ALL_B_14-25_mrc2 -- slowest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | lisbon_marathon | POR | 5 | 2021-25 | +0.0210 | [+0.0115, +0.0323] | +3.8 | 1.00 | [1, 4] | headline |
| 2 | nyc_marathon | USA | 11 | 2014-25 | +0.0143 | [+0.0125, +0.0151] | +2.6 | 1.00 | [2, 4] | headline |
| 3 | boston_marathon | USA | 11 | 2014-25 | +0.0139 | [+0.0125, +0.0160] | +2.5 | 1.00 | [1, 4] | headline |
| 4 | rome_marathon | ITA | 5 | 2021-25 | +0.0119 | [+0.0100, +0.0138] | +2.1 | 1.00 | [3, 6] | headline |
| 5 | madrid_marathon | ESP | 5 | 2021-25 | +0.0098 | [+0.0062, +0.0138] | +1.8 | 1.00 | [3, 7] | headline |
| 6 | milton_keynes_marathon | GBR | 11 | 2014-25 | +0.0079 | [+0.0032, +0.0148] | +1.4 | 1.00 | [2, 10] | headline |
| 7 | melbourne_marathon | AUS | 11 | 2014-25 | +0.0059 | [+0.0031, +0.0104] | +1.1 | 1.00 | [5, 11] | headline |
| 8 | edinburgh_marathon | GBR | 10 | 2014-25 | +0.0057 | [+0.0027, +0.0082] | +1.0 | 1.00 | [6, 11] | headline |
| 9 | oslo_marathon | NOR | 8 | 2016-25 | +0.0042 | [+0.0011, +0.0078] | +0.8 | 1.00 | [6, 12] | headline |
| 10 | vienna_marathon | AUT | 7 | 2018-25 | +0.0036 | [-0.0000, +0.0079] | +0.7 | 1.00 | [6, 13] | headline |
| 11 | belfast_marathon | GBR | 10 | 2014-25 | +0.0023 | [-0.0009, +0.0080] | +0.4 | 1.00 | [6, 14] | headline |
| 12 | copenhagen_marathon | DEN | 10 | 2014-25 | +0.0016 | [-0.0000, +0.0032] | +0.3 | 1.00 | [9, 13] | headline |
| 13 | london_marathon | GBR | 12 | 2014-25 | -0.0010 | [-0.0036, +0.0014] | -0.2 | 1.00 | [11, 16] | headline |
| 14 | munich_marathon | GER | 7 | 2018-25 | -0.0014 | [-0.0031, +0.0022] | -0.3 | 1.00 | [11, 15] | headline |
| 15 | cologne_marathon | GER | 6 | 2018-25 | -0.0024 | [-0.0044, -0.0001] | -0.4 | 1.00 | [12, 18] | headline |
| 16 | chester_marathon | GBR | 9 | 2016-25 | -0.0047 | [-0.0072, -0.0020] | -0.8 | 0.91 | [15, 22] | headline |
| 17 | prague_marathon | CZE | 10 | 2014-25 | -0.0053 | [-0.0098, +0.0000] | -0.9 | 0.66 | [14, 29] | headline |
| 18 | yorkshire_marathon | GBR | 6 | 2019-25 | -0.0057 | [-0.0084, -0.0026] | -1.0 | 0.76 | [15, 26] | headline |
| 19 | tokyo_marathon | JPN | 9 | 2015-25 | -0.0061 | [-0.0084, -0.0046] | -1.1 | 0.50 | [17, 26] | headline |
| 20 | hamburg_marathon | GER | 7 | 2018-25 | -0.0064 | [-0.0070, -0.0042] | -1.1 | 0.90 | [16, 21] | headline |
| 21 | stockholm_marathon | SWE | 11 | 2014-25 | -0.0070 | [-0.0095, -0.0049] | -1.3 | 0.34 | [17, 28] | tie |
| 23 | barcelona_marathon | ESP | 5 | 2021-25 | -0.0076 | [-0.0106, -0.0047] | -1.4 | 0.31 | [17, 30] | tie |
| 26 | malaga_marathon | ESP | 7 | 2018-25 | -0.0083 | [-0.0112, -0.0025] | -1.5 | 0.29 | [15, 31] | tie |
