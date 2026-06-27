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

slice: ALL_W_14-25_mrc2; min_n = 0; min_editions = 5; n_top = 20; model = full_nu8p00

### ALL_W_14-25_mrc2 -- fastest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | newport_marathon | GBR | 6 | 2018-25 | -0.0157 | [-0.0202, -0.0093] | -2.8 | 1.00 | [1, 7] | headline |
| 2 | dublin_marathon | IRL | 10 | 2014-25 | -0.0147 | [-0.0210, -0.0111] | -2.6 | 1.00 | [1, 6] | headline |
| 3 | amsterdam_marathon | NED | 9 | 2016-25 | -0.0130 | [-0.0153, -0.0092] | -2.3 | 1.00 | [1, 7] | headline |
| 4 | eindhoven_marathon | NED | 8 | 2017-25 | -0.0119 | [-0.0176, -0.0058] | -2.1 | 1.00 | [1, 15] | headline |
| 5 | seville_marathon | ESP | 7 | 2018-25 | -0.0115 | [-0.0160, -0.0068] | -2.1 | 1.00 | [1, 13] | headline |
| 6 | frankfurt_marathon | GER | 9 | 2015-25 | -0.0107 | [-0.0133, -0.0081] | -1.9 | 1.00 | [3, 9] | headline |
| 7 | malaga_marathon | ESP | 7 | 2018-25 | -0.0097 | [-0.0143, -0.0053] | -1.8 | 1.00 | [2, 16] | headline |
| 8 | paris_marathon | FRA | 11 | 2014-25 | -0.0092 | [-0.0114, -0.0046] | -1.7 | 0.98 | [4, 19] | headline |
| 9 | rotterdam_marathon | NED | 7 | 2018-25 | -0.0082 | [-0.0096, -0.0049] | -1.5 | 1.00 | [7, 19] | headline |
| 10 | manchester_marathon | GBR | 6 | 2019-25 | -0.0080 | [-0.0110, -0.0061] | -1.4 | 1.00 | [5, 16] | headline |
| 11 | cologne_marathon | GER | 6 | 2018-25 | -0.0076 | [-0.0111, -0.0019] | -1.4 | 0.86 | [5, 23] | headline |
| 12 | valencia_marathon | ESP | 9 | 2017-25 | -0.0072 | [-0.0101, -0.0045] | -1.3 | 1.00 | [6, 19] | headline |
| 13 | chester_marathon | GBR | 9 | 2016-25 | -0.0071 | [-0.0129, -0.0029] | -1.3 | 0.94 | [3, 21] | headline |
| 14 | barcelona_marathon | ESP | 5 | 2021-25 | -0.0070 | [-0.0114, -0.0021] | -1.3 | 0.89 | [6, 24] | headline |
| 15 | zurich_marathon | SUI | 5 | 2019-25 | -0.0067 | [-0.0122, -0.0013] | -1.2 | 0.82 | [3, 23] | headline |
| 16 | bostonUK_marathon | GBR | 5 | 2021-25 | -0.0055 | [-0.0123, +0.0014] | -1.0 | 0.66 | [4, 27] | headline |
| 17 | chicago_marathon | USA | 11 | 2014-25 | -0.0054 | [-0.0070, -0.0034] | -1.0 | 0.97 | [12, 21] | headline |
| 18 | berlin_marathon | GER | 11 | 2014-25 | -0.0053 | [-0.0067, -0.0030] | -1.0 | 0.87 | [13, 22] | headline |
| 19 | hannover_marathon | GER | 6 | 2018-25 | -0.0053 | [-0.0098, -0.0010] | -1.0 | 0.78 | [6, 24] | headline |
| 20 | stockholm_marathon | SWE | 11 | 2014-25 | -0.0051 | [-0.0103, +0.0010] | -0.9 | 0.65 | [7, 26] | headline |
| 21 | yorkshire_marathon | GBR | 6 | 2019-25 | -0.0050 | [-0.0101, -0.0010] | -0.9 | 0.67 | [6, 24] | headline |
| 22 | tokyo_marathon | JPN | 9 | 2015-25 | -0.0028 | [-0.0062, +0.0010] | -0.5 | 0.34 | [15, 25] | tie |
| 23 | hamburg_marathon | GER | 7 | 2018-25 | -0.0021 | [-0.0054, -0.0002] | -0.4 | 0.35 | [15, 25] | tie |

### ALL_W_14-25_mrc2 -- slowest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | madrid_marathon | ESP | 5 | 2021-25 | +0.0157 | [+0.0064, +0.0273] | +2.8 | 1.00 | [1, 5] | headline |
| 2 | nyc_marathon | USA | 11 | 2014-25 | +0.0136 | [+0.0118, +0.0150] | +2.4 | 1.00 | [1, 3] | headline |
| 3 | boston_marathon | USA | 11 | 2014-25 | +0.0065 | [+0.0041, +0.0100] | +1.2 | 1.00 | [3, 9] | headline |
| 4 | edinburgh_marathon | GBR | 10 | 2014-25 | +0.0060 | [+0.0015, +0.0088] | +1.1 | 1.00 | [3, 11] | headline |
| 5 | vienna_marathon | AUT | 7 | 2018-25 | +0.0057 | [+0.0019, +0.0139] | +1.0 | 1.00 | [2, 11] | headline |
| 6 | milton_keynes_marathon | GBR | 11 | 2014-25 | +0.0043 | [-0.0024, +0.0154] | +0.8 | 1.00 | [1, 15] | headline |
| 7 | copenhagen_marathon | DEN | 10 | 2014-25 | +0.0032 | [+0.0012, +0.0065] | +0.6 | 1.00 | [5, 12] | headline |
| 8 | oslo_marathon | NOR | 8 | 2016-25 | +0.0031 | [-0.0071, +0.0060] | +0.6 | 0.96 | [5, 23] | headline |
| 9 | munich_marathon | GER | 7 | 2018-25 | +0.0029 | [-0.0010, +0.0105] | +0.5 | 1.00 | [3, 13] | headline |
| 10 | melbourne_marathon | AUS | 11 | 2014-25 | +0.0029 | [-0.0010, +0.0100] | +0.5 | 1.00 | [3, 13] | headline |
| 11 | belfast_marathon | GBR | 10 | 2014-25 | +0.0020 | [-0.0046, +0.0116] | +0.4 | 0.99 | [3, 18] | headline |
| 12 | london_marathon | GBR | 11 | 2014-25 | +0.0017 | [+0.0005, +0.0036] | +0.3 | 1.00 | [8, 12] | headline |
| 13 | prague_marathon | CZE | 10 | 2014-25 | +0.0002 | [-0.0048, +0.0061] | +0.0 | 0.98 | [6, 20] | headline |
| 14 | hamburg_marathon | GER | 7 | 2018-25 | -0.0021 | [-0.0054, -0.0002] | -0.4 | 0.96 | [12, 22] | headline |
| 15 | tokyo_marathon | JPN | 9 | 2015-25 | -0.0028 | [-0.0062, +0.0010] | -0.5 | 0.94 | [12, 22] | headline |
| 16 | yorkshire_marathon | GBR | 6 | 2019-25 | -0.0050 | [-0.0101, -0.0010] | -0.9 | 0.59 | [13, 31] | headline |
| 17 | stockholm_marathon | SWE | 11 | 2014-25 | -0.0051 | [-0.0103, +0.0010] | -0.9 | 0.68 | [11, 30] | headline |
| 18 | hannover_marathon | GER | 6 | 2018-25 | -0.0053 | [-0.0098, -0.0010] | -1.0 | 0.48 | [13, 31] | tie |
| 19 | berlin_marathon | GER | 11 | 2014-25 | -0.0053 | [-0.0067, -0.0030] | -1.0 | 0.68 | [15, 24] | headline |
| 20 | chicago_marathon | USA | 11 | 2014-25 | -0.0054 | [-0.0070, -0.0034] | -1.0 | 0.58 | [16, 25] | headline |
| 21 | bostonUK_marathon | GBR | 5 | 2021-25 | -0.0055 | [-0.0123, +0.0014] | -1.0 | 0.54 | [10, 33] | headline |
| 22 | zurich_marathon | SUI | 5 | 2019-25 | -0.0067 | [-0.0122, -0.0013] | -1.2 | 0.35 | [14, 34] | tie |
| 23 | barcelona_marathon | ESP | 5 | 2021-25 | -0.0070 | [-0.0114, -0.0021] | -1.3 | 0.36 | [13, 31] | tie |
| 26 | cologne_marathon | GER | 6 | 2018-25 | -0.0076 | [-0.0111, -0.0019] | -1.4 | 0.37 | [14, 32] | tie |
