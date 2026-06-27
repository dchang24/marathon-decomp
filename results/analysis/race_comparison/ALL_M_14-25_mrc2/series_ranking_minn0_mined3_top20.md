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

slice: ALL_M_14-25_mrc2; min_n = 0; min_editions = 3; n_top = 20; model = full_nu8p00

### ALL_M_14-25_mrc2 -- fastest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | zurich_marathon | SUI | 5 | 2019-25 | -0.0182 | [-0.0206, -0.0144] | -3.2 | 1.00 | [1, 6] | headline |
| 2 | newport_marathon | GBR | 6 | 2018-25 | -0.0180 | [-0.0213, -0.0133] | -3.2 | 1.00 | [1, 8] | headline |
| 3 | seville_marathon | ESP | 7 | 2018-25 | -0.0170 | [-0.0189, -0.0147] | -3.0 | 1.00 | [1, 6] | headline |
| 4 | hannover_marathon | GER | 6 | 2018-25 | -0.0164 | [-0.0183, -0.0133] | -2.9 | 1.00 | [1, 9] | headline |
| 5 | eindhoven_marathon | NED | 8 | 2017-25 | -0.0144 | [-0.0161, -0.0095] | -2.6 | 1.00 | [3, 15] | headline |
| 6 | manchester_marathon | GBR | 6 | 2019-25 | -0.0142 | [-0.0159, -0.0125] | -2.5 | 1.00 | [4, 11] | headline |
| 7 | frankfurt_marathon | GER | 9 | 2015-25 | -0.0142 | [-0.0159, -0.0126] | -2.5 | 1.00 | [4, 11] | headline |
| 8 | valencia_marathon | ESP | 9 | 2017-25 | -0.0131 | [-0.0139, -0.0114] | -2.3 | 1.00 | [6, 13] | headline |
| 9 | dublin_marathon | IRL | 10 | 2014-25 | -0.0130 | [-0.0164, -0.0091] | -2.3 | 1.00 | [2, 16] | headline |
| 10 | barcelona_marathon | ESP | 5 | 2021-25 | -0.0125 | [-0.0153, -0.0093] | -2.2 | 1.00 | [4, 16] | headline |
| 11 | yorkshire_marathon | GBR | 6 | 2019-25 | -0.0124 | [-0.0145, -0.0079] | -2.2 | 0.99 | [6, 19] | headline |
| 12 | paris_marathon | FRA | 11 | 2014-25 | -0.0119 | [-0.0136, -0.0064] | -2.1 | 0.98 | [7, 20] | headline |
| 13 | amsterdam_marathon | NED | 9 | 2016-25 | -0.0107 | [-0.0129, -0.0090] | -1.9 | 1.00 | [8, 16] | headline |
| 14 | berlin_marathon | GER | 11 | 2014-25 | -0.0103 | [-0.0116, -0.0093] | -1.8 | 1.00 | [11, 16] | headline |
| 15 | bostonUK_marathon | GBR | 5 | 2021-25 | -0.0101 | [-0.0176, -0.0070] | -1.8 | 0.97 | [2, 21] | headline |
| 16 | san_sebastian_marathon | ESP | 3 | 2022-25 | -0.0092 | [-0.0139, -0.0015] | -1.6 | 0.78 | [8, 27] | headline |
| 17 | chester_marathon | GBR | 9 | 2016-25 | -0.0091 | [-0.0105, -0.0030] | -1.6 | 0.52 | [14, 26] | headline |
| 18 | hamburg_marathon | GER | 7 | 2018-25 | -0.0082 | [-0.0096, -0.0060] | -1.5 | 0.86 | [15, 22] | headline |
| 19 | rotterdam_marathon | NED | 7 | 2018-25 | -0.0076 | [-0.0095, -0.0056] | -1.4 | 0.75 | [16, 23] | headline |
| 20 | abingdon_marathon | GBR | 4 | 2022-25 | -0.0072 | [-0.0094, -0.0029] | -1.3 | 0.34 | [15, 27] | tie |
| 21 | malaga_marathon | ESP | 7 | 2018-25 | -0.0069 | [-0.0144, -0.0029] | -1.2 | 0.61 | [6, 27] | headline |
| 22 | chicago_marathon | USA | 11 | 2014-25 | -0.0068 | [-0.0081, -0.0050] | -1.2 | 0.60 | [17, 24] | headline |
| 24 | prague_marathon | CZE | 10 | 2014-25 | -0.0048 | [-0.0107, -0.0011] | -0.9 | 0.27 | [13, 27] | tie |

### ALL_M_14-25_mrc2 -- slowest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | venice_marathon | ITA | 3 | 2023-25 | +0.0271 | [+0.0234, +0.0292] | +4.9 | 1.00 | [1, 2] | headline |
| 2 | lisbon_marathon | POR | 5 | 2021-25 | +0.0208 | [+0.0112, +0.0318] | +3.8 | 1.00 | [1, 8] | headline |
| 3 | sydney_marathon | AUS | 4 | 2022-25 | +0.0196 | [+0.0160, +0.0229] | +3.6 | 1.00 | [2, 5] | headline |
| 4 | nyc_marathon | USA | 11 | 2014-25 | +0.0175 | [+0.0160, +0.0188] | +3.2 | 1.00 | [3, 6] | headline |
| 5 | melbourne_marathon | AUS | 11 | 2014-25 | +0.0159 | [+0.0100, +0.0194] | +2.9 | 1.00 | [3, 9] | headline |
| 6 | milton_keynes_marathon | GBR | 11 | 2014-25 | +0.0152 | [+0.0095, +0.0189] | +2.7 | 1.00 | [3, 10] | headline |
| 7 | cape_town_marathon | RSA | 3 | 2022-24 | +0.0149 | [+0.0025, +0.0264] | +2.7 | 1.00 | [1, 15] | headline |
| 8 | brighton_marathon | GBR | 4 | 2022-25 | +0.0104 | [+0.0067, +0.0145] | +1.9 | 1.00 | [5, 11] | headline |
| 9 | boston_marathon | USA | 11 | 2014-25 | +0.0097 | [+0.0090, +0.0130] | +1.8 | 1.00 | [7, 10] | headline |
| 10 | belfast_marathon | GBR | 10 | 2014-25 | +0.0097 | [+0.0033, +0.0148] | +1.8 | 1.00 | [5, 14] | headline |
| 11 | oslo_marathon | NOR | 8 | 2016-25 | +0.0066 | [+0.0033, +0.0092] | +1.2 | 1.00 | [10, 15] | headline |
| 12 | milan_marathon | ITA | 4 | 2022-25 | +0.0064 | [+0.0039, +0.0083] | +1.2 | 1.00 | [9, 15] | headline |
| 13 | madrid_marathon | ESP | 5 | 2021-25 | +0.0054 | [+0.0020, +0.0104] | +1.0 | 1.00 | [8, 16] | headline |
| 14 | rome_marathon | ITA | 5 | 2021-25 | +0.0053 | [+0.0036, +0.0083] | +1.0 | 1.00 | [10, 15] | headline |
| 15 | edinburgh_marathon | GBR | 10 | 2014-25 | +0.0022 | [-0.0008, +0.0059] | +0.4 | 1.00 | [12, 18] | headline |
| 16 | vienna_marathon | AUT | 7 | 2018-25 | +0.0012 | [-0.0016, +0.0057] | +0.2 | 1.00 | [13, 19] | headline |
| 17 | copenhagen_marathon | DEN | 10 | 2014-25 | +0.0003 | [-0.0023, +0.0022] | +0.1 | 1.00 | [15, 19] | headline |
| 18 | tokyo_marathon | JPN | 9 | 2015-25 | -0.0002 | [-0.0029, +0.0011] | -0.0 | 0.98 | [16, 20] | headline |
| 19 | london_marathon | GBR | 12 | 2014-25 | -0.0044 | [-0.0051, -0.0009] | -0.8 | 0.50 | [18, 24] | headline |
| 20 | cologne_marathon | GER | 6 | 2018-25 | -0.0045 | [-0.0068, -0.0013] | -0.8 | 0.33 | [18, 26] | tie |
| 21 | stockholm_marathon | SWE | 11 | 2014-25 | -0.0046 | [-0.0070, +0.0009] | -0.8 | 0.50 | [17, 27] | headline |
