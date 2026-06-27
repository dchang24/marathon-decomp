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

slice: ALL_B_14-25_mrc2; min_n = 0; min_editions = 3; n_top = 20; model = full_nu8p00

### ALL_B_14-25_mrc2 -- fastest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | newport_marathon | GBR | 6 | 2018-25 | -0.0173 | [-0.0208, -0.0130] | -3.1 | 1.00 | [1, 6] | headline |
| 2 | zurich_marathon | SUI | 5 | 2019-25 | -0.0169 | [-0.0192, -0.0127] | -3.0 | 1.00 | [1, 7] | headline |
| 3 | seville_marathon | ESP | 7 | 2018-25 | -0.0161 | [-0.0171, -0.0144] | -2.9 | 1.00 | [1, 5] | headline |
| 4 | hannover_marathon | GER | 6 | 2018-25 | -0.0143 | [-0.0162, -0.0112] | -2.6 | 1.00 | [2, 12] | headline |
| 5 | dublin_marathon | IRL | 10 | 2014-25 | -0.0142 | [-0.0174, -0.0106] | -2.5 | 1.00 | [1, 14] | headline |
| 6 | eindhoven_marathon | NED | 8 | 2017-25 | -0.0139 | [-0.0167, -0.0103] | -2.5 | 1.00 | [1, 15] | headline |
| 7 | manchester_marathon | GBR | 6 | 2019-25 | -0.0138 | [-0.0152, -0.0123] | -2.5 | 1.00 | [3, 9] | headline |
| 8 | frankfurt_marathon | GER | 9 | 2015-25 | -0.0128 | [-0.0144, -0.0111] | -2.3 | 1.00 | [4, 13] | headline |
| 9 | barcelona_marathon | ESP | 5 | 2021-25 | -0.0123 | [-0.0148, -0.0096] | -2.2 | 1.00 | [4, 17] | headline |
| 10 | valencia_marathon | ESP | 9 | 2017-25 | -0.0118 | [-0.0127, -0.0099] | -2.1 | 1.00 | [8, 16] | headline |
| 11 | yorkshire_marathon | GBR | 6 | 2019-25 | -0.0114 | [-0.0129, -0.0078] | -2.0 | 0.98 | [7, 19] | headline |
| 12 | bostonUK_marathon | GBR | 5 | 2021-25 | -0.0114 | [-0.0153, -0.0080] | -2.0 | 1.00 | [3, 18] | headline |
| 13 | paris_marathon | FRA | 11 | 2014-25 | -0.0113 | [-0.0133, -0.0066] | -2.0 | 0.98 | [6, 20] | headline |
| 14 | amsterdam_marathon | NED | 9 | 2016-25 | -0.0108 | [-0.0129, -0.0095] | -1.9 | 1.00 | [8, 16] | headline |
| 15 | san_sebastian_marathon | ESP | 3 | 2022-25 | -0.0106 | [-0.0122, -0.0027] | -1.9 | 0.87 | [8, 26] | headline |
| 16 | chester_marathon | GBR | 9 | 2016-25 | -0.0103 | [-0.0125, -0.0049] | -1.8 | 0.92 | [9, 23] | headline |
| 17 | malaga_marathon | ESP | 7 | 2018-25 | -0.0092 | [-0.0157, -0.0046] | -1.7 | 0.81 | [4, 23] | headline |
| 18 | berlin_marathon | GER | 11 | 2014-25 | -0.0092 | [-0.0099, -0.0074] | -1.7 | 0.97 | [15, 21] | headline |
| 19 | rotterdam_marathon | NED | 7 | 2018-25 | -0.0081 | [-0.0094, -0.0063] | -1.4 | 0.89 | [15, 21] | headline |
| 20 | abingdon_marathon | GBR | 4 | 2022-25 | -0.0076 | [-0.0097, -0.0038] | -1.4 | 0.45 | [15, 25] | tie |
| 21 | hamburg_marathon | GER | 7 | 2018-25 | -0.0074 | [-0.0087, -0.0055] | -1.3 | 0.61 | [17, 22] | headline |
| 22 | cologne_marathon | GER | 6 | 2018-25 | -0.0064 | [-0.0082, -0.0041] | -1.1 | 0.40 | [17, 24] | tie |

### ALL_B_14-25_mrc2 -- slowest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | venice_marathon | ITA | 3 | 2023-25 | +0.0263 | [+0.0233, +0.0287] | +4.8 | 1.00 | [1, 2] | headline |
| 2 | HK_marathon | HKG | 4 | 2018-25 | +0.0180 | [+0.0147, +0.0247] | +3.3 | 1.00 | [2, 6] | headline |
| 3 | sydney_marathon | AUS | 4 | 2022-25 | +0.0179 | [+0.0145, +0.0207] | +3.3 | 1.00 | [2, 5] | headline |
| 4 | lisbon_marathon | POR | 5 | 2021-25 | +0.0170 | [+0.0089, +0.0282] | +3.1 | 1.00 | [1, 9] | headline |
| 5 | milton_keynes_marathon | GBR | 11 | 2014-25 | +0.0168 | [+0.0111, +0.0186] | +3.0 | 1.00 | [3, 7] | headline |
| 6 | nyc_marathon | USA | 11 | 2014-25 | +0.0143 | [+0.0139, +0.0159] | +2.6 | 1.00 | [4, 6] | headline |
| 7 | brighton_marathon | GBR | 4 | 2022-25 | +0.0103 | [+0.0076, +0.0126] | +1.9 | 1.00 | [6, 11] | headline |
| 8 | melbourne_marathon | AUS | 11 | 2014-25 | +0.0085 | [+0.0043, +0.0123] | +1.5 | 0.99 | [6, 15] | headline |
| 9 | boston_marathon | USA | 11 | 2014-25 | +0.0074 | [+0.0062, +0.0098] | +1.3 | 1.00 | [7, 13] | headline |
| 10 | milan_marathon | ITA | 4 | 2022-25 | +0.0071 | [+0.0049, +0.0090] | +1.3 | 1.00 | [8, 14] | headline |
| 11 | rome_marathon | ITA | 5 | 2021-25 | +0.0067 | [+0.0049, +0.0090] | +1.2 | 1.00 | [8, 14] | headline |
| 12 | madrid_marathon | ESP | 5 | 2021-25 | +0.0059 | [+0.0029, +0.0107] | +1.1 | 1.00 | [7, 15] | headline |
| 13 | belfast_marathon | GBR | 10 | 2014-25 | +0.0059 | [+0.0023, +0.0112] | +1.1 | 1.00 | [7, 16] | headline |
| 14 | oslo_marathon | NOR | 8 | 2016-25 | +0.0043 | [+0.0008, +0.0076] | +0.8 | 0.99 | [11, 18] | headline |
| 15 | edinburgh_marathon | GBR | 10 | 2014-25 | +0.0026 | [-0.0003, +0.0047] | +0.5 | 1.00 | [14, 19] | headline |
| 16 | vienna_marathon | AUT | 7 | 2018-25 | +0.0021 | [-0.0004, +0.0064] | +0.4 | 0.99 | [12, 18] | headline |
| 17 | cape_town_marathon | RSA | 3 | 2022-24 | +0.0015 | [-0.0068, +0.0111] | +0.3 | 0.84 | [7, 26] | headline |
| 18 | copenhagen_marathon | DEN | 10 | 2014-25 | +0.0011 | [-0.0006, +0.0023] | +0.2 | 1.00 | [15, 19] | headline |
| 19 | london_marathon | GBR | 12 | 2014-25 | -0.0011 | [-0.0031, -0.0002] | -0.2 | 0.66 | [18, 22] | headline |
| 20 | tokyo_marathon | JPN | 9 | 2015-25 | -0.0016 | [-0.0029, +0.0002] | -0.3 | 0.77 | [18, 22] | headline |
| 21 | stockholm_marathon | SWE | 11 | 2014-25 | -0.0025 | [-0.0062, +0.0005] | -0.4 | 0.40 | [18, 26] | tie |
| 22 | prague_marathon | CZE | 10 | 2014-25 | -0.0034 | [-0.0071, +0.0022] | -0.6 | 0.32 | [16, 27] | tie |
