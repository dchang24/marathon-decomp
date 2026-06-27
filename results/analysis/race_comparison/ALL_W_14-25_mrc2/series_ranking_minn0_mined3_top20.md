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

slice: ALL_W_14-25_mrc2; min_n = 0; min_editions = 3; n_top = 20; model = full_nu8p00

### ALL_W_14-25_mrc2 -- fastest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | newport_marathon | GBR | 6 | 2018-25 | -0.0173 | [-0.0211, -0.0097] | -3.1 | 1.00 | [1, 13] | headline |
| 2 | dublin_marathon | IRL | 10 | 2014-25 | -0.0165 | [-0.0200, -0.0106] | -2.9 | 1.00 | [1, 9] | headline |
| 3 | malaga_marathon | ESP | 7 | 2018-25 | -0.0144 | [-0.0176, -0.0088] | -2.6 | 1.00 | [1, 15] | headline |
| 4 | frankfurt_marathon | GER | 9 | 2015-25 | -0.0130 | [-0.0146, -0.0093] | -2.3 | 1.00 | [2, 12] | headline |
| 5 | seville_marathon | ESP | 7 | 2018-25 | -0.0128 | [-0.0163, -0.0099] | -2.3 | 1.00 | [1, 12] | headline |
| 6 | manchester_marathon | GBR | 6 | 2019-25 | -0.0119 | [-0.0143, -0.0094] | -2.1 | 1.00 | [3, 13] | headline |
| 7 | cologne_marathon | GER | 6 | 2018-25 | -0.0115 | [-0.0136, -0.0068] | -2.1 | 1.00 | [3, 18] | headline |
| 8 | barcelona_marathon | ESP | 5 | 2021-25 | -0.0114 | [-0.0148, -0.0043] | -2.0 | 0.98 | [3, 20] | headline |
| 9 | chester_marathon | GBR | 9 | 2016-25 | -0.0111 | [-0.0154, -0.0043] | -2.0 | 0.96 | [2, 21] | headline |
| 10 | amsterdam_marathon | NED | 9 | 2016-25 | -0.0106 | [-0.0141, -0.0081] | -1.9 | 1.00 | [4, 16] | headline |
| 11 | eindhoven_marathon | NED | 8 | 2017-25 | -0.0103 | [-0.0173, -0.0060] | -1.8 | 1.00 | [1, 19] | headline |
| 12 | zurich_marathon | SUI | 5 | 2019-25 | -0.0094 | [-0.0164, -0.0012] | -1.7 | 0.85 | [2, 25] | headline |
| 13 | rotterdam_marathon | NED | 7 | 2018-25 | -0.0090 | [-0.0117, -0.0060] | -1.6 | 1.00 | [7, 19] | headline |
| 14 | abingdon_marathon | GBR | 4 | 2022-25 | -0.0084 | [-0.0122, -0.0017] | -1.5 | 0.87 | [6, 24] | headline |
| 15 | yorkshire_marathon | GBR | 6 | 2019-25 | -0.0081 | [-0.0140, -0.0040] | -1.5 | 0.94 | [3, 22] | headline |
| 16 | paris_marathon | FRA | 11 | 2014-25 | -0.0080 | [-0.0116, -0.0034] | -1.4 | 0.91 | [8, 22] | headline |
| 17 | valencia_marathon | ESP | 9 | 2017-25 | -0.0078 | [-0.0114, -0.0037] | -1.4 | 0.96 | [8, 21] | headline |
| 18 | bostonUK_marathon | GBR | 5 | 2021-25 | -0.0071 | [-0.0139, -0.0020] | -1.3 | 0.83 | [3, 24] | headline |
| 19 | hannover_marathon | GER | 6 | 2018-25 | -0.0062 | [-0.0109, -0.0016] | -1.1 | 0.87 | [8, 24] | headline |
| 20 | hamburg_marathon | GER | 7 | 2018-25 | -0.0047 | [-0.0074, -0.0012] | -0.8 | 0.46 | [16, 26] | tie |
| 21 | berlin_marathon | GER | 11 | 2014-25 | -0.0046 | [-0.0071, -0.0031] | -0.8 | 0.61 | [15, 22] | headline |
| 22 | cape_town_marathon | RSA | 3 | 2022-24 | -0.0045 | [-0.0206, +0.0050] | -0.8 | 0.51 | [1, 33] | headline |

### ALL_W_14-25_mrc2 -- slowest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | venice_marathon | ITA | 3 | 2023-25 | +0.0231 | [+0.0184, +0.0289] | +4.2 | 1.00 | [1, 2] | headline |
| 2 | sydney_marathon | AUS | 4 | 2022-25 | +0.0163 | [+0.0120, +0.0215] | +3.0 | 1.00 | [1, 5] | headline |
| 3 | nyc_marathon | USA | 11 | 2014-25 | +0.0145 | [+0.0132, +0.0165] | +2.6 | 1.00 | [2, 5] | headline |
| 4 | madrid_marathon | ESP | 5 | 2021-25 | +0.0141 | [+0.0044, +0.0244] | +2.6 | 1.00 | [1, 12] | headline |
| 5 | milton_keynes_marathon | GBR | 11 | 2014-25 | +0.0107 | [+0.0031, +0.0186] | +1.9 | 1.00 | [2, 14] | headline |
| 6 | brighton_marathon | GBR | 4 | 2022-25 | +0.0107 | [+0.0069, +0.0158] | +1.9 | 1.00 | [3, 9] | headline |
| 7 | rome_marathon | ITA | 4 | 2021-25 | +0.0099 | [+0.0066, +0.0131] | +1.8 | 1.00 | [4, 10] | headline |
| 8 | milan_marathon | ITA | 4 | 2022-25 | +0.0087 | [+0.0050, +0.0123] | +1.6 | 1.00 | [5, 11] | headline |
| 9 | prague_marathon | CZE | 10 | 2014-25 | +0.0042 | [-0.0015, +0.0143] | +0.8 | 0.97 | [4, 21] | headline |
| 10 | vienna_marathon | AUT | 7 | 2018-25 | +0.0041 | [-0.0013, +0.0138] | +0.7 | 0.98 | [4, 20] | headline |
| 11 | boston_marathon | USA | 11 | 2014-25 | +0.0040 | [+0.0024, +0.0074] | +0.7 | 1.00 | [9, 14] | headline |
| 12 | edinburgh_marathon | GBR | 10 | 2014-25 | +0.0039 | [-0.0017, +0.0070] | +0.7 | 0.99 | [8, 20] | headline |
| 13 | copenhagen_marathon | DEN | 10 | 2014-25 | +0.0038 | [+0.0010, +0.0072] | +0.7 | 1.00 | [9, 17] | headline |
| 14 | stockholm_marathon | SWE | 11 | 2014-25 | +0.0027 | [-0.0031, +0.0086] | +0.5 | 0.95 | [8, 22] | headline |
| 15 | london_marathon | GBR | 11 | 2014-25 | +0.0003 | [-0.0021, +0.0019] | +0.1 | 0.94 | [14, 21] | headline |
| 16 | oslo_marathon | NOR | 8 | 2016-25 | +0.0002 | [-0.0064, +0.0051] | +0.0 | 0.72 | [11, 28] | headline |
| 17 | belfast_marathon | GBR | 10 | 2014-25 | +0.0002 | [-0.0050, +0.0118] | +0.0 | 0.78 | [6, 26] | headline |
| 18 | melbourne_marathon | AUS | 11 | 2014-25 | -0.0006 | [-0.0035, +0.0082] | -0.1 | 0.86 | [7, 23] | headline |
| 19 | tokyo_marathon | JPN | 9 | 2015-25 | -0.0007 | [-0.0026, +0.0011] | -0.1 | 0.84 | [15, 22] | headline |
| 20 | chicago_marathon | USA | 11 | 2014-25 | -0.0013 | [-0.0030, +0.0004] | -0.2 | 0.72 | [17, 23] | headline |
| 21 | munich_marathon | GER | 7 | 2018-25 | -0.0028 | [-0.0053, +0.0066] | -0.5 | 0.71 | [9, 24] | headline |
| 22 | cape_town_marathon | RSA | 3 | 2022-24 | -0.0045 | [-0.0206, +0.0050] | -0.8 | 0.31 | [11, 43] | tie |
