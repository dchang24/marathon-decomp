# Top-N slowest / fastest race editions by v_j

Estimand: v_j under the beta=0 (bundling) gauge = race effect RELATIVE
TO CONTEMPORANEOUS RACES; it bundles course + that day's weather +
field conditions. Selection is by bootstrap rank stability P(top-N)
(headline P>=0.5, tie P>=0.25), which demotes small-field races whose
extreme point v_j is noise (winner's curse). min@3:00 =
180*(exp(v_j)-1) = minutes vs the average race for a 3:00:00 runner.

slice: Po10_B_14-25_mrc2; min_n = 0; n_top = 20; model = full_nu8p00

### Po10_B_14-25_mrc2 -- slowest 20 (headline P>=0.5, tie P>=0.25)

| rank | race | country | n_j | v_j | v 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|
| 1 | milton_keynes_marathon 2018 | GBR | 533 | +0.1097 | [+0.1013, +0.1167] | +20.9 | 1.00 | [1, 2] | headline |
| 2 | stockholm_marathon 2024 | SWE | 32 | +0.0995 | [+0.0838, +0.1195] | +18.8 | 1.00 | [1, 3] | headline |
| 3 | nyc_marathon 2022 | USA | 582 | +0.0860 | [+0.0807, +0.0909] | +16.2 | 1.00 | [2, 4] | headline |
| 4 | vienna_marathon 2018 | AUT | 32 | +0.0670 | [+0.0385, +0.0910] | +12.5 | 0.99 | [3, 15] | headline |
| 5 | london_marathon 2018 | GBR | 6947 | +0.0654 | [+0.0631, +0.0665] | +12.2 | 1.00 | [4, 6] | headline |
| 6 | athens_marathon 2025 | GRE | 63 | +0.0592 | [+0.0395, +0.0804] | +11.0 | 1.00 | [4, 14] | headline |
| 7 | boston_marathon 2017 | USA | 275 | +0.0557 | [+0.0487, +0.0629] | +10.3 | 1.00 | [5, 9] | headline |
| 8 | berlin_marathon 2025 | GER | 1048 | +0.0473 | [+0.0450, +0.0527] | +8.7 | 1.00 | [7, 12] | headline |
| 9 | boston_marathon 2024 | USA | 652 | +0.0473 | [+0.0443, +0.0509] | +8.7 | 1.00 | [7, 13] | headline |
| 10 | venice_marathon 2024 | ITA | 33 | +0.0436 | [+0.0227, +0.0721] | +8.0 | 0.82 | [4, 36] | headline |
| 11 | stockholm_marathon 2022 | SWE | 39 | +0.0430 | [+0.0208, +0.0636] | +7.9 | 0.89 | [5, 39] | headline |
| 12 | chicago_marathon 2021 | USA | 38 | +0.0425 | [+0.0126, +0.0706] | +7.8 | 0.68 | [4, 64] | headline |
| 13 | copenhagen_marathon 2016 | DEN | 55 | +0.0390 | [+0.0240, +0.0592] | +7.2 | 0.79 | [7, 32] | headline |
| 14 | boston_marathon 2016 | USA | 209 | +0.0385 | [+0.0314, +0.0464] | +7.1 | 0.93 | [10, 23] | headline |
| 15 | copenhagen_marathon 2018 | DEN | 79 | +0.0376 | [+0.0217, +0.0519] | +6.9 | 0.71 | [8, 36] | headline |
| 16 | vienna_marathon 2023 | AUT | 64 | +0.0372 | [+0.0232, +0.0517] | +6.8 | 0.69 | [7, 34] | headline |
| 17 | boston_marathon 2018 | USA | 320 | +0.0342 | [+0.0297, +0.0420] | +6.3 | 0.78 | [11, 26] | headline |
| 18 | boston_marathon 2019 | USA | 465 | +0.0336 | [+0.0290, +0.0386] | +6.1 | 0.68 | [13, 26] | headline |
| 19 | nyc_marathon 2023 | USA | 549 | +0.0327 | [+0.0280, +0.0373] | +6.0 | 0.52 | [13, 28] | headline |
| 20 | stockholm_marathon 2023 | SWE | 43 | +0.0325 | [+0.0201, +0.0469] | +6.0 | 0.61 | [10, 40] | headline |
| 21 | athens_marathon 2024 | GRE | 83 | +0.0325 | [+0.0201, +0.0532] | +6.0 | 0.54 | [8, 41] | headline |
| 22 | milton_keynes_marathon 2023 | GBR | 440 | +0.0312 | [+0.0252, +0.0357] | +5.7 | 0.29 | [15, 33] | tie |
| 23 | chicago_marathon 2017 | USA | 287 | +0.0298 | [+0.0218, +0.0360] | +5.4 | 0.28 | [16, 37] | tie |
| 27 | prague_marathon 2024 | CZE | 33 | +0.0260 | [+0.0071, +0.0446] | +4.7 | 0.26 | [12, 78] | tie |
| 28 | munich_marathon 2023 | GER | 35 | +0.0254 | [+0.0013, +0.0455] | +4.6 | 0.29 | [11, 104] | tie |

### Po10_B_14-25_mrc2 -- fastest 20 (headline P>=0.5, tie P>=0.25)

| rank | race | country | n_j | v_j | v 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|
| 1 | hamburg_marathon 2019 | GER | 31 | -0.0422 | [-0.0552, -0.0238] | -7.4 | 0.97 | [1, 21] | headline |
| 2 | amsterdam_marathon 2025 | NED | 222 | -0.0350 | [-0.0396, -0.0264] | -6.2 | 0.98 | [2, 15] | headline |
| 3 | rotterdam_marathon 2021 | NED | 23 | -0.0329 | [-0.0481, -0.0150] | -5.8 | 0.85 | [1, 71] | headline |
| 4 | seville_marathon 2019 | ESP | 89 | -0.0306 | [-0.0406, -0.0215] | -5.4 | 0.90 | [2, 27] | headline |
| 5 | malaga_marathon 2019 | ESP | 53 | -0.0298 | [-0.0418, -0.0170] | -5.3 | 0.74 | [1, 56] | headline |
| 6 | amsterdam_marathon 2023 | NED | 347 | -0.0293 | [-0.0349, -0.0255] | -5.2 | 0.98 | [4, 18] | headline |
| 7 | munich_marathon 2025 | GER | 29 | -0.0291 | [-0.0641, +0.0020] | -5.2 | 0.61 | [1, 165] | headline |
| 8 | london_marathon 2015 | GBR | 7230 | -0.0286 | [-0.0306, -0.0259] | -5.1 | 0.99 | [6, 17] | headline |
| 9 | berlin_marathon 2014 | GER | 517 | -0.0272 | [-0.0315, -0.0207] | -4.8 | 0.76 | [5, 34] | headline |
| 10 | tokyo_marathon 2015 | JPN | 47 | -0.0265 | [-0.0342, -0.0146] | -4.7 | 0.66 | [4, 75] | headline |
| 11 | paris_marathon 2025 | FRA | 99 | -0.0264 | [-0.0371, -0.0163] | -4.7 | 0.61 | [2, 56] | headline |
| 12 | copenhagen_marathon 2025 | DEN | 105 | -0.0260 | [-0.0326, -0.0144] | -4.6 | 0.51 | [5, 71] | headline |
| 13 | barcelona_marathon 2025 | ESP | 302 | -0.0258 | [-0.0299, -0.0187] | -4.6 | 0.68 | [6, 42] | headline |
| 14 | dublin_marathon 2019 | IRL | 173 | -0.0258 | [-0.0341, -0.0185] | -4.6 | 0.67 | [3, 49] | headline |
| 15 | san_sebastian_marathon 2025 | ESP | 49 | -0.0257 | [-0.0446, -0.0080] | -4.6 | 0.56 | [2, 111] | headline |
| 16 | paris_marathon 2019 | FRA | 223 | -0.0255 | [-0.0342, -0.0187] | -4.5 | 0.64 | [4, 44] | headline |
| 17 | chicago_marathon 2014 | USA | 126 | -0.0255 | [-0.0330, -0.0144] | -4.5 | 0.50 | [4, 77] | headline |
| 18 | barcelona_marathon 2021 | ESP | 73 | -0.0252 | [-0.0365, -0.0145] | -4.5 | 0.69 | [2, 73] | headline |
| 19 | valencia_marathon 2023 | ESP | 812 | -0.0249 | [-0.0281, -0.0215] | -4.4 | 0.63 | [10, 30] | headline |
| 20 | boston_marathon 2015 | USA | 191 | -0.0236 | [-0.0303, -0.0174] | -4.2 | 0.36 | [7, 53] | tie |
| 22 | seville_marathon 2020 | ESP | 114 | -0.0221 | [-0.0335, -0.0139] | -3.9 | 0.31 | [4, 78] | tie |
| 26 | newport_marathon 2023 | GBR | 479 | -0.0217 | [-0.0273, -0.0168] | -3.9 | 0.27 | [11, 56] | tie |
| 27 | amsterdam_marathon 2019 | NED | 278 | -0.0216 | [-0.0286, -0.0162] | -3.9 | 0.40 | [8, 59] | tie |
| 28 | frankfurt_marathon 2025 | GER | 160 | -0.0213 | [-0.0299, -0.0117] | -3.8 | 0.27 | [8, 88] | tie |
| 29 | melbourne_marathon 2025 | AUS | 20 | -0.0212 | [-0.0504, +0.0033] | -3.8 | 0.49 | [1, 171] | tie |
| 30 | vienna_marathon 2025 | AUT | 89 | -0.0210 | [-0.0334, -0.0068] | -3.7 | 0.27 | [4, 122] | tie |
| 31 | newport_marathon 2019 | GBR | 500 | -0.0207 | [-0.0279, -0.0141] | -3.7 | 0.26 | [11, 75] | tie |
