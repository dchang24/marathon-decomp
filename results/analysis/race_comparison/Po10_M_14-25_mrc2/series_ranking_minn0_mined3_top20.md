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

slice: Po10_M_14-25_mrc2; min_n = 0; min_editions = 3; n_top = 20; model = full_nu8p00

### Po10_M_14-25_mrc2 -- fastest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | barcelona_marathon | ESP | 4 | 2021-25 | -0.0215 | [-0.0257, -0.0133] | -3.8 | 1.00 | [1, 8] | headline |
| 2 | seville_marathon | ESP | 7 | 2018-25 | -0.0186 | [-0.0226, -0.0128] | -3.3 | 1.00 | [1, 10] | headline |
| 3 | newport_marathon | GBR | 6 | 2018-25 | -0.0176 | [-0.0202, -0.0123] | -3.1 | 1.00 | [1, 9] | headline |
| 4 | valencia_marathon | ESP | 8 | 2017-25 | -0.0175 | [-0.0197, -0.0120] | -3.1 | 1.00 | [1, 10] | headline |
| 5 | dublin_marathon | IRL | 10 | 2014-25 | -0.0165 | [-0.0181, -0.0109] | -2.9 | 1.00 | [2, 11] | headline |
| 6 | frankfurt_marathon | GER | 9 | 2015-25 | -0.0156 | [-0.0193, -0.0101] | -2.8 | 1.00 | [1, 11] | headline |
| 7 | amsterdam_marathon | NED | 9 | 2016-25 | -0.0151 | [-0.0184, -0.0063] | -2.7 | 1.00 | [2, 15] | headline |
| 8 | manchester_marathon | GBR | 6 | 2019-25 | -0.0145 | [-0.0152, -0.0118] | -2.6 | 1.00 | [4, 10] | headline |
| 9 | paris_marathon | FRA | 11 | 2014-25 | -0.0129 | [-0.0176, -0.0057] | -2.3 | 1.00 | [2, 15] | headline |
| 10 | yorkshire_marathon | GBR | 6 | 2019-25 | -0.0117 | [-0.0140, -0.0065] | -2.1 | 1.00 | [7, 16] | headline |
| 11 | hamburg_marathon | GER | 5 | 2018-25 | -0.0105 | [-0.0185, +0.0097] | -1.9 | 0.79 | [1, 26] | headline |
| 12 | berlin_marathon | GER | 11 | 2014-25 | -0.0105 | [-0.0120, -0.0080] | -1.9 | 1.00 | [9, 14] | headline |
| 13 | bostonUK_marathon | GBR | 5 | 2021-25 | -0.0104 | [-0.0189, -0.0052] | -1.9 | 1.00 | [2, 16] | headline |
| 14 | chicago_marathon | USA | 11 | 2014-25 | -0.0100 | [-0.0124, -0.0039] | -1.8 | 1.00 | [8, 16] | headline |
| 15 | chester_marathon | GBR | 9 | 2016-25 | -0.0082 | [-0.0101, -0.0012] | -1.5 | 1.00 | [10, 19] | headline |
| 16 | abingdon_marathon | GBR | 4 | 2022-25 | -0.0075 | [-0.0098, -0.0027] | -1.3 | 1.00 | [12, 19] | headline |
| 17 | malaga_marathon | ESP | 5 | 2018-25 | -0.0045 | [-0.0112, +0.0034] | -0.8 | 0.85 | [10, 24] | headline |
| 18 | london_marathon | GBR | 11 | 2014-25 | -0.0029 | [-0.0041, -0.0006] | -0.5 | 0.94 | [16, 21] | headline |
| 19 | prague_marathon | CZE | 4 | 2022-25 | -0.0026 | [-0.0127, +0.0201] | -0.5 | 0.46 | [7, 30] | tie |
| 21 | copenhagen_marathon | DEN | 9 | 2015-25 | +0.0017 | [-0.0080, +0.0082] | +0.3 | 0.50 | [13, 27] | headline |
| 22 | rotterdam_marathon | NED | 5 | 2018-24 | +0.0018 | [-0.0075, +0.0083] | +0.3 | 0.43 | [14, 27] | tie |
| 23 | vienna_marathon | AUT | 6 | 2018-25 | +0.0023 | [-0.0055, +0.0194] | +0.4 | 0.28 | [16, 31] | tie |

### Po10_M_14-25_mrc2 -- slowest 20 series (headline P>=0.5, tie P>=0.25)

| rank | series | country | k | years | median v | m 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | boston_marathon | USA | 11 | 2014-25 | +0.0300 | [+0.0142, +0.0338] | +5.5 | 1.00 | [1, 4] | headline |
| 2 | nyc_marathon | USA | 11 | 2014-25 | +0.0214 | [+0.0180, +0.0254] | +3.9 | 1.00 | [1, 4] | headline |
| 3 | munich_marathon | GER | 3 | 2023-25 | +0.0146 | [-0.0084, +0.0350] | +2.6 | 0.97 | [1, 21] | headline |
| 4 | milton_keynes_marathon | GBR | 11 | 2014-25 | +0.0121 | [+0.0075, +0.0191] | +2.2 | 1.00 | [2, 7] | headline |
| 5 | brighton_marathon | GBR | 4 | 2022-25 | +0.0108 | [+0.0086, +0.0149] | +2.0 | 1.00 | [3, 7] | headline |
| 6 | stockholm_marathon | SWE | 8 | 2014-24 | +0.0105 | [+0.0027, +0.0256] | +1.9 | 1.00 | [2, 10] | headline |
| 7 | belfast_marathon | GBR | 10 | 2014-25 | +0.0042 | [-0.0027, +0.0100] | +0.8 | 1.00 | [6, 15] | headline |
| 8 | rome_marathon | ITA | 4 | 2022-25 | +0.0037 | [-0.0051, +0.0127] | +0.7 | 1.00 | [4, 17] | headline |
| 9 | tokyo_marathon | JPN | 8 | 2015-25 | +0.0024 | [-0.0022, +0.0061] | +0.4 | 1.00 | [7, 15] | headline |
| 10 | vienna_marathon | AUT | 6 | 2018-25 | +0.0023 | [-0.0055, +0.0194] | +0.4 | 1.00 | [2, 17] | headline |
| 11 | rotterdam_marathon | NED | 5 | 2018-24 | +0.0018 | [-0.0075, +0.0083] | +0.3 | 0.99 | [6, 19] | headline |
| 12 | copenhagen_marathon | DEN | 9 | 2015-25 | +0.0017 | [-0.0080, +0.0082] | +0.3 | 0.99 | [6, 20] | headline |
| 13 | edinburgh_marathon | GBR | 10 | 2014-25 | +0.0013 | [-0.0022, +0.0050] | +0.2 | 1.00 | [8, 16] | headline |
| 14 | prague_marathon | CZE | 4 | 2022-25 | -0.0026 | [-0.0127, +0.0201] | -0.5 | 0.94 | [3, 26] | headline |
| 15 | london_marathon | GBR | 11 | 2014-25 | -0.0029 | [-0.0041, -0.0006] | -0.5 | 1.00 | [12, 17] | headline |
| 16 | malaga_marathon | ESP | 5 | 2018-25 | -0.0045 | [-0.0112, +0.0034] | -0.8 | 0.94 | [9, 23] | headline |
| 17 | abingdon_marathon | GBR | 4 | 2022-25 | -0.0075 | [-0.0098, -0.0027] | -1.3 | 0.92 | [14, 21] | headline |
| 18 | chester_marathon | GBR | 9 | 2016-25 | -0.0082 | [-0.0101, -0.0012] | -1.5 | 0.92 | [14, 23] | headline |
| 19 | chicago_marathon | USA | 11 | 2014-25 | -0.0100 | [-0.0124, -0.0039] | -1.8 | 0.74 | [17, 25] | headline |
| 22 | hamburg_marathon | GER | 5 | 2018-25 | -0.0105 | [-0.0185, +0.0097] | -1.9 | 0.65 | [7, 32] | headline |
| 23 | yorkshire_marathon | GBR | 6 | 2019-25 | -0.0117 | [-0.0140, -0.0065] | -2.1 | 0.31 | [17, 26] | tie |
