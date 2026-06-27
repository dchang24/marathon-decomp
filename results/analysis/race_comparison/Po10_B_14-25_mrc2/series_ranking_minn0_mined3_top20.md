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

slice: Po10_B_14-25_mrc2; min_n = 0; min_editions = 3; n_top = 20; model = full_nu8p00

### Po10_B_14-25_mrc2 -- fastest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | barcelona_marathon | ESP | 5 | 2021-25 | -0.0189 | [-0.0241, -0.0105] | -3.4 | 1.00 | [1, 10] | headline |
| 2 | dublin_marathon | IRL | 10 | 2014-25 | -0.0174 | [-0.0191, -0.0120] | -3.1 | 1.00 | [1, 9] | headline |
| 3 | seville_marathon | ESP | 7 | 2018-25 | -0.0174 | [-0.0205, -0.0115] | -3.1 | 1.00 | [1, 8] | headline |
| 4 | newport_marathon | GBR | 6 | 2018-25 | -0.0171 | [-0.0199, -0.0122] | -3.1 | 1.00 | [1, 9] | headline |
| 5 | frankfurt_marathon | GER | 9 | 2015-25 | -0.0169 | [-0.0200, -0.0107] | -3.0 | 1.00 | [1, 11] | headline |
| 6 | valencia_marathon | ESP | 8 | 2017-25 | -0.0162 | [-0.0169, -0.0126] | -2.9 | 1.00 | [2, 8] | headline |
| 7 | manchester_marathon | GBR | 6 | 2019-25 | -0.0136 | [-0.0145, -0.0116] | -2.4 | 1.00 | [5, 10] | headline |
| 8 | paris_marathon | FRA | 11 | 2014-25 | -0.0128 | [-0.0174, -0.0084] | -2.3 | 1.00 | [2, 13] | headline |
| 9 | bostonUK_marathon | GBR | 5 | 2021-25 | -0.0116 | [-0.0154, -0.0065] | -2.1 | 1.00 | [4, 16] | headline |
| 10 | yorkshire_marathon | GBR | 6 | 2019-25 | -0.0111 | [-0.0129, -0.0064] | -2.0 | 1.00 | [7, 16] | headline |
| 11 | amsterdam_marathon | NED | 9 | 2016-25 | -0.0101 | [-0.0156, -0.0073] | -1.8 | 1.00 | [5, 15] | headline |
| 12 | chester_marathon | GBR | 9 | 2016-25 | -0.0099 | [-0.0117, -0.0049] | -1.8 | 1.00 | [10, 17] | headline |
| 13 | berlin_marathon | GER | 11 | 2014-25 | -0.0082 | [-0.0101, -0.0061] | -1.5 | 1.00 | [9, 17] | headline |
| 14 | malaga_marathon | ESP | 7 | 2018-25 | -0.0081 | [-0.0142, -0.0004] | -1.5 | 0.98 | [6, 20] | headline |
| 15 | chicago_marathon | USA | 11 | 2014-25 | -0.0080 | [-0.0091, -0.0022] | -1.4 | 1.00 | [11, 19] | headline |
| 16 | abingdon_marathon | GBR | 4 | 2022-25 | -0.0076 | [-0.0098, -0.0028] | -1.4 | 1.00 | [10, 20] | headline |
| 17 | cologne_marathon | GER | 4 | 2019-25 | -0.0040 | [-0.0164, +0.0085] | -0.7 | 0.73 | [3, 27] | headline |
| 18 | rotterdam_marathon | NED | 6 | 2018-24 | -0.0037 | [-0.0121, +0.0030] | -0.7 | 0.86 | [10, 25] | headline |
| 19 | hamburg_marathon | GER | 6 | 2018-25 | -0.0033 | [-0.0135, +0.0078] | -0.6 | 0.71 | [6, 26] | headline |
| 20 | prague_marathon | CZE | 4 | 2022-25 | -0.0016 | [-0.0085, +0.0139] | -0.3 | 0.38 | [12, 29] | tie |
| 21 | london_marathon | GBR | 11 | 2014-25 | -0.0013 | [-0.0026, +0.0003] | -0.2 | 0.62 | [18, 22] | headline |

### Po10_B_14-25_mrc2 -- slowest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | boston_marathon | USA | 11 | 2014-25 | +0.0266 | [+0.0151, +0.0346] | +4.9 | 1.00 | [1, 5] | headline |
| 2 | nyc_marathon | USA | 11 | 2014-25 | +0.0207 | [+0.0179, +0.0239] | +3.8 | 1.00 | [1, 4] | headline |
| 3 | munich_marathon | GER | 4 | 2022-25 | +0.0127 | [+0.0007, +0.0245] | +2.3 | 1.00 | [1, 13] | headline |
| 4 | stockholm_marathon | SWE | 8 | 2014-24 | +0.0121 | [+0.0075, +0.0263] | +2.2 | 1.00 | [1, 8] | headline |
| 5 | brighton_marathon | GBR | 4 | 2022-25 | +0.0109 | [+0.0089, +0.0141] | +2.0 | 1.00 | [4, 8] | headline |
| 6 | vienna_marathon | AUT | 6 | 2018-25 | +0.0108 | [+0.0001, +0.0216] | +2.0 | 1.00 | [2, 12] | headline |
| 7 | milton_keynes_marathon | GBR | 11 | 2014-25 | +0.0094 | [+0.0054, +0.0166] | +1.7 | 1.00 | [3, 9] | headline |
| 8 | copenhagen_marathon | DEN | 9 | 2015-25 | +0.0059 | [-0.0051, +0.0097] | +1.1 | 1.00 | [6, 17] | headline |
| 9 | rome_marathon | ITA | 4 | 2022-25 | +0.0025 | [-0.0064, +0.0099] | +0.5 | 1.00 | [6, 19] | headline |
| 10 | belfast_marathon | GBR | 10 | 2014-25 | +0.0024 | [-0.0029, +0.0082] | +0.4 | 1.00 | [8, 15] | headline |
| 11 | edinburgh_marathon | GBR | 10 | 2014-25 | +0.0019 | [-0.0006, +0.0058] | +0.3 | 1.00 | [8, 14] | headline |
| 12 | tokyo_marathon | JPN | 8 | 2015-25 | +0.0011 | [-0.0021, +0.0042] | +0.2 | 1.00 | [8, 15] | headline |
| 13 | london_marathon | GBR | 11 | 2014-25 | -0.0013 | [-0.0026, +0.0003] | -0.2 | 1.00 | [12, 16] | headline |
| 14 | prague_marathon | CZE | 4 | 2022-25 | -0.0016 | [-0.0085, +0.0139] | -0.3 | 0.96 | [5, 22] | headline |
| 15 | hamburg_marathon | GER | 6 | 2018-25 | -0.0033 | [-0.0135, +0.0078] | -0.6 | 0.81 | [8, 28] | headline |
| 16 | rotterdam_marathon | NED | 6 | 2018-24 | -0.0037 | [-0.0121, +0.0030] | -0.7 | 0.80 | [9, 24] | headline |
| 17 | cologne_marathon | GER | 4 | 2019-25 | -0.0040 | [-0.0164, +0.0085] | -0.7 | 0.74 | [7, 31] | headline |
| 18 | abingdon_marathon | GBR | 4 | 2022-25 | -0.0076 | [-0.0098, -0.0028] | -1.4 | 0.84 | [14, 24] | headline |
| 19 | chicago_marathon | USA | 11 | 2014-25 | -0.0080 | [-0.0091, -0.0022] | -1.4 | 0.91 | [15, 23] | headline |
| 20 | malaga_marathon | ESP | 7 | 2018-25 | -0.0081 | [-0.0142, -0.0004] | -1.5 | 0.68 | [14, 28] | headline |
| 21 | berlin_marathon | GER | 11 | 2014-25 | -0.0082 | [-0.0101, -0.0061] | -1.5 | 0.45 | [17, 25] | tie |
| 22 | chester_marathon | GBR | 9 | 2016-25 | -0.0099 | [-0.0117, -0.0049] | -1.8 | 0.41 | [17, 24] | tie |
