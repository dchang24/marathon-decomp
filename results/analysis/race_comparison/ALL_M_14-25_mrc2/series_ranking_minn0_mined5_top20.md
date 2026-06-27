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

slice: ALL_M_14-25_mrc2; min_n = 0; min_editions = 5; n_top = 20; model = full_nu8p00

### ALL_M_14-25_mrc2 -- fastest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | frankfurt_marathon | GER | 9 | 2015-25 | -0.0166 | [-0.0183, -0.0149] | -3.0 | 1.00 | [1, 3] | headline |
| 2 | newport_marathon | GBR | 6 | 2018-25 | -0.0158 | [-0.0204, -0.0130] | -2.8 | 1.00 | [1, 6] | headline |
| 3 | hannover_marathon | GER | 6 | 2018-25 | -0.0142 | [-0.0160, -0.0113] | -2.6 | 1.00 | [2, 11] | headline |
| 4 | zurich_marathon | SUI | 5 | 2019-25 | -0.0135 | [-0.0170, -0.0096] | -2.4 | 1.00 | [2, 13] | headline |
| 5 | amsterdam_marathon | NED | 9 | 2016-25 | -0.0132 | [-0.0146, -0.0111] | -2.4 | 1.00 | [3, 10] | headline |
| 6 | seville_marathon | ESP | 7 | 2018-25 | -0.0131 | [-0.0168, -0.0109] | -2.4 | 1.00 | [2, 9] | headline |
| 7 | eindhoven_marathon | NED | 8 | 2017-25 | -0.0129 | [-0.0160, -0.0090] | -2.3 | 1.00 | [2, 14] | headline |
| 8 | dublin_marathon | IRL | 10 | 2014-25 | -0.0122 | [-0.0160, -0.0084] | -2.2 | 1.00 | [2, 17] | headline |
| 9 | valencia_marathon | ESP | 9 | 2017-25 | -0.0116 | [-0.0135, -0.0099] | -2.1 | 1.00 | [5, 12] | headline |
| 10 | paris_marathon | FRA | 11 | 2014-25 | -0.0112 | [-0.0138, -0.0092] | -2.0 | 1.00 | [5, 15] | headline |
| 11 | manchester_marathon | GBR | 6 | 2019-25 | -0.0105 | [-0.0120, -0.0088] | -1.9 | 1.00 | [5, 15] | headline |
| 12 | berlin_marathon | GER | 11 | 2014-25 | -0.0102 | [-0.0110, -0.0084] | -1.8 | 1.00 | [10, 15] | headline |
| 13 | prague_marathon | CZE | 10 | 2014-25 | -0.0087 | [-0.0158, -0.0032] | -1.6 | 0.90 | [2, 23] | headline |
| 14 | bostonUK_marathon | GBR | 5 | 2021-25 | -0.0086 | [-0.0157, -0.0032] | -1.6 | 0.88 | [2, 24] | headline |
| 15 | chicago_marathon | USA | 11 | 2014-25 | -0.0081 | [-0.0097, -0.0071] | -1.5 | 1.00 | [12, 18] | headline |
| 16 | barcelona_marathon | ESP | 5 | 2021-25 | -0.0080 | [-0.0106, -0.0049] | -1.4 | 0.96 | [11, 21] | headline |
| 17 | tokyo_marathon | JPN | 9 | 2015-25 | -0.0070 | [-0.0091, -0.0030] | -1.3 | 0.71 | [14, 24] | headline |
| 18 | rotterdam_marathon | NED | 7 | 2018-25 | -0.0069 | [-0.0090, -0.0050] | -1.2 | 0.87 | [14, 22] | headline |
| 19 | yorkshire_marathon | GBR | 6 | 2019-25 | -0.0061 | [-0.0098, -0.0021] | -1.1 | 0.61 | [13, 24] | headline |
| 20 | stockholm_marathon | SWE | 11 | 2014-25 | -0.0060 | [-0.0085, -0.0035] | -1.1 | 0.66 | [14, 23] | headline |
| 21 | malaga_marathon | ESP | 7 | 2018-25 | -0.0060 | [-0.0114, -0.0003] | -1.1 | 0.62 | [9, 26] | headline |
| 22 | hamburg_marathon | GER | 7 | 2018-25 | -0.0059 | [-0.0076, -0.0037] | -1.1 | 0.64 | [17, 23] | headline |

### ALL_M_14-25_mrc2 -- slowest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | lisbon_marathon | POR | 5 | 2021-25 | +0.0237 | [+0.0152, +0.0360] | +4.3 | 1.00 | [1, 3] | headline |
| 2 | boston_marathon | USA | 11 | 2014-25 | +0.0161 | [+0.0150, +0.0192] | +2.9 | 1.00 | [1, 3] | headline |
| 3 | nyc_marathon | USA | 11 | 2014-25 | +0.0147 | [+0.0120, +0.0157] | +2.6 | 1.00 | [3, 6] | headline |
| 4 | milton_keynes_marathon | GBR | 11 | 2014-25 | +0.0134 | [+0.0048, +0.0178] | +2.4 | 1.00 | [2, 9] | headline |
| 5 | rome_marathon | ITA | 5 | 2021-25 | +0.0115 | [+0.0091, +0.0138] | +2.1 | 1.00 | [4, 7] | headline |
| 6 | melbourne_marathon | AUS | 11 | 2014-25 | +0.0108 | [+0.0067, +0.0163] | +1.9 | 1.00 | [3, 8] | headline |
| 7 | madrid_marathon | ESP | 5 | 2021-25 | +0.0098 | [+0.0060, +0.0135] | +1.8 | 1.00 | [3, 9] | headline |
| 8 | belfast_marathon | GBR | 10 | 2014-25 | +0.0064 | [-0.0007, +0.0105] | +1.2 | 1.00 | [6, 13] | headline |
| 9 | oslo_marathon | NOR | 8 | 2016-25 | +0.0056 | [+0.0034, +0.0091] | +1.0 | 1.00 | [7, 11] | headline |
| 10 | edinburgh_marathon | GBR | 10 | 2014-25 | +0.0055 | [+0.0011, +0.0082] | +1.0 | 1.00 | [8, 12] | headline |
| 11 | vienna_marathon | AUT | 7 | 2018-25 | +0.0027 | [-0.0014, +0.0070] | +0.5 | 1.00 | [8, 14] | headline |
| 12 | copenhagen_marathon | DEN | 10 | 2014-25 | +0.0010 | [-0.0012, +0.0029] | +0.2 | 1.00 | [10, 14] | headline |
| 13 | cologne_marathon | GER | 6 | 2018-25 | -0.0007 | [-0.0032, +0.0023] | -0.1 | 1.00 | [11, 15] | headline |
| 14 | munich_marathon | GER | 7 | 2018-25 | -0.0020 | [-0.0040, +0.0021] | -0.4 | 1.00 | [11, 17] | headline |
| 15 | chester_marathon | GBR | 9 | 2016-25 | -0.0023 | [-0.0071, +0.0004] | -0.4 | 0.96 | [12, 21] | headline |
| 16 | london_marathon | GBR | 12 | 2014-25 | -0.0037 | [-0.0049, -0.0005] | -0.7 | 1.00 | [13, 18] | headline |
| 17 | hamburg_marathon | GER | 7 | 2018-25 | -0.0059 | [-0.0076, -0.0037] | -1.1 | 0.75 | [16, 22] | headline |
| 18 | malaga_marathon | ESP | 7 | 2018-25 | -0.0060 | [-0.0114, -0.0003] | -1.1 | 0.49 | [13, 30] | tie |
| 19 | stockholm_marathon | SWE | 11 | 2014-25 | -0.0060 | [-0.0085, -0.0035] | -1.1 | 0.64 | [16, 25] | headline |
| 20 | yorkshire_marathon | GBR | 6 | 2019-25 | -0.0061 | [-0.0098, -0.0021] | -1.1 | 0.59 | [15, 26] | headline |
| 21 | rotterdam_marathon | NED | 7 | 2018-25 | -0.0069 | [-0.0090, -0.0050] | -1.2 | 0.40 | [17, 25] | tie |
| 22 | tokyo_marathon | JPN | 9 | 2015-25 | -0.0070 | [-0.0091, -0.0030] | -1.3 | 0.58 | [15, 25] | headline |
