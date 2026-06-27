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

slice: WA_M_14-25_mrc2; min_n = 0; min_editions = 3; n_top = 20; model = full_nu8p00

### WA_M_14-25_mrc2 -- fastest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | hannover_marathon | GER | 3 | 2022-25 | -0.0191 | [-0.0239, -0.0120] | -3.4 | 1.00 | [1, 5] | headline |
| 2 | valencia_marathon | ESP | 9 | 2017-25 | -0.0134 | [-0.0158, -0.0121] | -2.4 | 1.00 | [1, 6] | headline |
| 3 | berlin_marathon | GER | 11 | 2014-25 | -0.0131 | [-0.0145, -0.0094] | -2.3 | 1.00 | [2, 6] | headline |
| 4 | zurich_marathon | SUI | 3 | 2023-25 | -0.0123 | [-0.0222, -0.0064] | -2.2 | 1.00 | [1, 9] | headline |
| 5 | seville_marathon | ESP | 7 | 2018-25 | -0.0123 | [-0.0190, -0.0100] | -2.2 | 1.00 | [1, 7] | headline |
| 6 | tokyo_marathon | JPN | 9 | 2015-25 | -0.0083 | [-0.0133, -0.0043] | -1.5 | 1.00 | [3, 12] | headline |
| 7 | eindhoven_marathon | NED | 5 | 2019-25 | -0.0077 | [-0.0134, -0.0028] | -1.4 | 1.00 | [3, 14] | headline |
| 8 | barcelona_marathon | ESP | 5 | 2021-25 | -0.0065 | [-0.0130, -0.0021] | -1.2 | 1.00 | [5, 14] | headline |
| 9 | frankfurt_marathon | GER | 9 | 2015-25 | -0.0059 | [-0.0099, -0.0015] | -1.1 | 1.00 | [6, 14] | headline |
| 10 | paris_marathon | FRA | 11 | 2014-25 | -0.0050 | [-0.0104, -0.0016] | -0.9 | 1.00 | [6, 14] | headline |
| 11 | melbourne_marathon | AUS | 4 | 2019-25 | -0.0048 | [-0.0140, +0.0095] | -0.9 | 0.92 | [4, 24] | headline |
| 12 | rotterdam_marathon | NED | 7 | 2018-25 | -0.0047 | [-0.0070, +0.0016] | -0.8 | 1.00 | [8, 18] | headline |
| 13 | amsterdam_marathon | NED | 9 | 2016-25 | -0.0041 | [-0.0103, -0.0014] | -0.7 | 1.00 | [6, 14] | headline |
| 14 | dublin_marathon | IRL | 9 | 2015-25 | -0.0037 | [-0.0086, +0.0037] | -0.7 | 0.97 | [7, 21] | headline |
| 15 | hamburg_marathon | GER | 6 | 2018-25 | -0.0009 | [-0.0036, +0.0046] | -0.2 | 0.94 | [12, 22] | headline |
| 16 | london_marathon | GBR | 12 | 2014-25 | +0.0005 | [-0.0013, +0.0031] | +0.1 | 0.96 | [14, 21] | headline |
| 17 | rome_marathon | ITA | 3 | 2023-25 | +0.0009 | [-0.0067, +0.0098] | +0.2 | 0.68 | [9, 25] | headline |
| 18 | chicago_marathon | USA | 11 | 2014-25 | +0.0022 | [-0.0027, +0.0048] | +0.4 | 0.90 | [13, 22] | headline |
| 19 | copenhagen_marathon | DEN | 5 | 2018-25 | +0.0024 | [-0.0017, +0.0106] | +0.4 | 0.56 | [13, 25] | headline |
| 20 | malaga_marathon | ESP | 3 | 2023-25 | +0.0031 | [-0.0029, +0.0102] | +0.6 | 0.57 | [13, 25] | headline |
| 21 | prague_marathon | CZE | 4 | 2019-25 | +0.0033 | [-0.0043, +0.0113] | +0.6 | 0.49 | [11, 25] | tie |
| 22 | milan_marathon | ITA | 4 | 2022-25 | +0.0049 | [-0.0005, +0.0110] | +0.9 | 0.38 | [14, 26] | tie |
| 23 | vienna_marathon | AUT | 5 | 2019-25 | +0.0054 | [-0.0002, +0.0100] | +1.0 | 0.34 | [16, 25] | tie |
| 24 | manchester_marathon | GBR | 6 | 2019-25 | +0.0072 | [-0.0017, +0.0108] | +1.3 | 0.29 | [16, 25] | tie |

### WA_M_14-25_mrc2 -- slowest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | nyc_marathon | USA | 11 | 2014-25 | +0.0223 | [+0.0185, +0.0262] | +4.1 | 1.00 | [1, 2] | headline |
| 2 | stockholm_marathon | SWE | 5 | 2021-25 | +0.0207 | [+0.0154, +0.0267] | +3.8 | 1.00 | [1, 2] | headline |
| 3 | sydney_marathon | AUS | 3 | 2023-25 | +0.0127 | [+0.0061, +0.0185] | +2.3 | 1.00 | [3, 6] | headline |
| 4 | boston_marathon | USA | 11 | 2014-25 | +0.0093 | [+0.0070, +0.0140] | +1.7 | 1.00 | [3, 7] | headline |
| 5 | manchester_marathon | GBR | 6 | 2019-25 | +0.0072 | [-0.0017, +0.0108] | +1.3 | 0.99 | [4, 13] | headline |
| 6 | vienna_marathon | AUT | 5 | 2019-25 | +0.0054 | [-0.0002, +0.0100] | +1.0 | 1.00 | [4, 13] | headline |
| 7 | milan_marathon | ITA | 4 | 2022-25 | +0.0049 | [-0.0005, +0.0110] | +0.9 | 1.00 | [3, 15] | headline |
| 8 | prague_marathon | CZE | 4 | 2019-25 | +0.0033 | [-0.0043, +0.0113] | +0.6 | 0.99 | [4, 18] | headline |
| 9 | malaga_marathon | ESP | 3 | 2023-25 | +0.0031 | [-0.0029, +0.0102] | +0.6 | 1.00 | [4, 16] | headline |
| 10 | copenhagen_marathon | DEN | 5 | 2018-25 | +0.0024 | [-0.0017, +0.0106] | +0.4 | 1.00 | [4, 16] | headline |
| 11 | chicago_marathon | USA | 11 | 2014-25 | +0.0022 | [-0.0027, +0.0048] | +0.4 | 1.00 | [7, 16] | headline |
| 12 | rome_marathon | ITA | 3 | 2023-25 | +0.0009 | [-0.0067, +0.0098] | +0.2 | 0.98 | [4, 20] | headline |
| 13 | london_marathon | GBR | 12 | 2014-25 | +0.0005 | [-0.0013, +0.0031] | +0.1 | 1.00 | [8, 15] | headline |
| 14 | hamburg_marathon | GER | 6 | 2018-25 | -0.0009 | [-0.0036, +0.0046] | -0.2 | 1.00 | [7, 17] | headline |
| 15 | dublin_marathon | IRL | 9 | 2015-25 | -0.0037 | [-0.0086, +0.0037] | -0.7 | 0.94 | [8, 22] | headline |
| 16 | amsterdam_marathon | NED | 9 | 2016-25 | -0.0041 | [-0.0103, -0.0014] | -0.7 | 0.74 | [15, 23] | headline |
| 17 | rotterdam_marathon | NED | 7 | 2018-25 | -0.0047 | [-0.0070, +0.0016] | -0.8 | 0.95 | [11, 21] | headline |
| 18 | melbourne_marathon | AUS | 4 | 2019-25 | -0.0048 | [-0.0140, +0.0095] | -0.9 | 0.69 | [5, 25] | headline |
| 19 | paris_marathon | FRA | 11 | 2014-25 | -0.0050 | [-0.0104, -0.0016] | -0.9 | 0.66 | [15, 23] | headline |
| 20 | frankfurt_marathon | GER | 9 | 2015-25 | -0.0059 | [-0.0099, -0.0015] | -1.1 | 0.71 | [15, 23] | headline |
| 21 | barcelona_marathon | ESP | 5 | 2021-25 | -0.0065 | [-0.0130, -0.0021] | -1.2 | 0.47 | [15, 24] | tie |
| 22 | eindhoven_marathon | NED | 5 | 2019-25 | -0.0077 | [-0.0134, -0.0028] | -1.4 | 0.53 | [15, 26] | headline |
| 23 | tokyo_marathon | JPN | 9 | 2015-25 | -0.0083 | [-0.0133, -0.0043] | -1.5 | 0.32 | [17, 26] | tie |
