# Top-N slowest / fastest race editions by v_j

Estimand: v_j under the beta=0 (bundling) gauge = race effect RELATIVE
TO CONTEMPORANEOUS RACES; it bundles course + that day's weather +
field conditions. Selection is by bootstrap rank stability P(top-N)
(headline P>=0.5, tie P>=0.25), which demotes small-field races whose
extreme point v_j is noise (winner's curse). min@3:00 =
180*(exp(v_j)-1) = minutes vs the average race for a 3:00:00 runner.

slice: Po10_W_14-25_mrc2; min_n = 0; n_top = 20; model = full_nu8p00

### Po10_W_14-25_mrc2 -- slowest 20 (headline P>=0.5, tie P>=0.25)

| rank | race | country | n_j | v_j | v 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|
| 1 | milton_keynes_marathon 2018 | GBR | 136 | +0.1133 | [+0.1015, +0.1248] | +21.6 | 1.00 | [1, 1] | headline |
| 2 | nyc_marathon 2022 | USA | 214 | +0.0791 | [+0.0705, +0.0849] | +14.8 | 1.00 | [2, 3] | headline |
| 3 | london_marathon 2018 | GBR | 2526 | +0.0636 | [+0.0600, +0.0663] | +11.8 | 1.00 | [3, 5] | headline |
| 4 | athens_marathon 2024 | GRE | 28 | +0.0528 | [+0.0161, +0.0818] | +9.8 | 0.94 | [2, 39] | headline |
| 5 | boston_marathon 2017 | USA | 85 | +0.0504 | [+0.0387, +0.0589] | +9.3 | 1.00 | [4, 11] | headline |
| 6 | athens_marathon 2025 | GRE | 23 | +0.0504 | [+0.0220, +0.0755] | +9.3 | 0.95 | [3, 28] | headline |
| 7 | berlin_marathon 2025 | GER | 393 | +0.0499 | [+0.0446, +0.0575] | +9.2 | 1.00 | [4, 9] | headline |
| 8 | boston_marathon 2024 | USA | 249 | +0.0438 | [+0.0381, +0.0483] | +8.1 | 1.00 | [5, 13] | headline |
| 9 | paris_marathon 2018 | FRA | 42 | +0.0423 | [+0.0290, +0.0560] | +7.8 | 0.99 | [4, 20] | headline |
| 10 | boston_marathon 2019 | USA | 162 | +0.0412 | [+0.0340, +0.0517] | +7.6 | 1.00 | [6, 16] | headline |
| 11 | copenhagen_marathon 2018 | DEN | 24 | +0.0412 | [+0.0188, +0.0645] | +7.6 | 0.90 | [4, 33] | headline |
| 12 | boston_marathon 2016 | USA | 78 | +0.0383 | [+0.0270, +0.0519] | +7.0 | 0.98 | [5, 19] | headline |
| 13 | nyc_marathon 2023 | USA | 191 | +0.0368 | [+0.0303, +0.0443] | +6.7 | 0.98 | [8, 19] | headline |
| 14 | chicago_marathon 2017 | USA | 90 | +0.0321 | [+0.0235, +0.0444] | +5.9 | 0.93 | [9, 24] | headline |
| 15 | nyc_marathon 2015 | USA | 102 | +0.0309 | [+0.0208, +0.0407] | +5.7 | 0.80 | [10, 30] | headline |
| 16 | boston_marathon 2018 | USA | 111 | +0.0292 | [+0.0180, +0.0415] | +5.3 | 0.66 | [10, 34] | headline |
| 17 | copenhagen_marathon 2023 | DEN | 32 | +0.0282 | [-0.0015, +0.0643] | +5.2 | 0.52 | [3, 92] | headline |
| 18 | belfast_marathon 2017 | GBR | 26 | +0.0272 | [+0.0035, +0.0491] | +5.0 | 0.55 | [6, 73] | headline |
| 19 | vienna_marathon 2023 | AUT | 26 | +0.0267 | [+0.0068, +0.0453] | +4.9 | 0.54 | [9, 58] | headline |
| 20 | milton_keynes_marathon 2023 | GBR | 143 | +0.0265 | [+0.0173, +0.0361] | +4.8 | 0.45 | [14, 38] | tie |
| 21 | edinburgh_marathon 2017 | GBR | 197 | +0.0248 | [+0.0162, +0.0363] | +4.5 | 0.32 | [13, 40] | tie |
| 22 | berlin_marathon 2021 | GER | 142 | +0.0245 | [+0.0152, +0.0342] | +4.5 | 0.29 | [14, 39] | tie |
| 23 | belfast_marathon 2018 | GBR | 47 | +0.0236 | [+0.0069, +0.0407] | +4.3 | 0.33 | [11, 64] | tie |
| 24 | nyc_marathon 2025 | USA | 246 | +0.0229 | [+0.0175, +0.0309] | +4.2 | 0.34 | [17, 35] | tie |
| 39 | copenhagen_marathon 2019 | DEN | 20 | +0.0158 | [-0.0234, +0.0589] | +2.9 | 0.31 | [4, 184] | tie |

### Po10_W_14-25_mrc2 -- fastest 20 (headline P>=0.5, tie P>=0.25)

| rank | race | country | n_j | v_j | v 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|
| 1 | belfast_marathon 2014 | GBR | 21 | -0.0374 | [-0.0629, -0.0115] | -6.6 | 0.85 | [1, 66] | headline |
| 2 | seville_marathon 2019 | ESP | 22 | -0.0336 | [-0.0524, -0.0082] | -5.9 | 0.86 | [1, 78] | headline |
| 3 | barcelona_marathon 2025 | ESP | 80 | -0.0300 | [-0.0399, -0.0204] | -5.3 | 0.95 | [1, 27] | headline |
| 4 | amsterdam_marathon 2025 | NED | 70 | -0.0295 | [-0.0399, -0.0179] | -5.2 | 0.86 | [1, 34] | headline |
| 5 | paris_marathon 2025 | FRA | 33 | -0.0294 | [-0.0485, -0.0083] | -5.2 | 0.82 | [1, 81] | headline |
| 6 | amsterdam_marathon 2019 | NED | 62 | -0.0263 | [-0.0432, -0.0157] | -4.7 | 0.77 | [1, 43] | headline |
| 7 | paris_marathon 2019 | FRA | 50 | -0.0255 | [-0.0484, -0.0077] | -4.5 | 0.66 | [1, 83] | headline |
| 8 | frankfurt_marathon 2025 | GER | 26 | -0.0252 | [-0.0341, -0.0117] | -4.5 | 0.64 | [2, 67] | headline |
| 9 | frankfurt_marathon 2016 | GER | 31 | -0.0244 | [-0.0459, +0.0029] | -4.3 | 0.61 | [2, 131] | headline |
| 10 | chester_marathon 2016 | GBR | 167 | -0.0243 | [-0.0332, -0.0151] | -4.3 | 0.65 | [3, 49] | headline |
| 11 | newport_marathon 2025 | GBR | 149 | -0.0242 | [-0.0331, -0.0168] | -4.3 | 0.72 | [4, 38] | headline |
| 12 | london_marathon 2015 | GBR | 2350 | -0.0240 | [-0.0265, -0.0206] | -4.3 | 0.79 | [9, 29] | headline |
| 13 | dublin_marathon 2019 | IRL | 56 | -0.0238 | [-0.0347, -0.0095] | -4.2 | 0.57 | [3, 74] | headline |
| 14 | dublin_marathon 2018 | IRL | 58 | -0.0237 | [-0.0355, -0.0132] | -4.2 | 0.59 | [3, 56] | headline |
| 15 | edinburgh_marathon 2014 | GBR | 160 | -0.0236 | [-0.0351, -0.0148] | -4.2 | 0.57 | [3, 53] | headline |
| 16 | amsterdam_marathon 2023 | NED | 122 | -0.0220 | [-0.0290, -0.0159] | -3.9 | 0.55 | [8, 45] | headline |
| 17 | paris_marathon 2015 | FRA | 91 | -0.0216 | [-0.0317, -0.0114] | -3.8 | 0.42 | [4, 69] | tie |
| 18 | chester_marathon 2025 | GBR | 208 | -0.0212 | [-0.0272, -0.0114] | -3.8 | 0.31 | [11, 66] | tie |
| 19 | valencia_marathon 2023 | ESP | 222 | -0.0211 | [-0.0270, -0.0166] | -3.8 | 0.31 | [10, 43] | tie |
| 20 | barcelona_marathon 2024 | ESP | 42 | -0.0210 | [-0.0362, -0.0080] | -3.7 | 0.47 | [4, 82] | tie |
| 21 | dublin_marathon 2015 | IRL | 74 | -0.0210 | [-0.0347, -0.0093] | -3.7 | 0.48 | [3, 78] | tie |
| 22 | boston_marathon 2015 | USA | 75 | -0.0201 | [-0.0273, -0.0103] | -3.6 | 0.30 | [9, 70] | tie |
| 23 | frankfurt_marathon 2019 | GER | 39 | -0.0196 | [-0.0353, -0.0066] | -3.5 | 0.39 | [3, 90] | tie |
| 25 | newport_marathon 2019 | GBR | 148 | -0.0192 | [-0.0305, -0.0070] | -3.4 | 0.31 | [6, 88] | tie |
| 30 | chicago_marathon 2014 | USA | 36 | -0.0183 | [-0.0335, -0.0023] | -3.3 | 0.28 | [5, 104] | tie |
| 31 | copenhagen_marathon 2025 | DEN | 20 | -0.0175 | [-0.0299, +0.0015] | -3.1 | 0.27 | [7, 125] | tie |
| 46 | seville_marathon 2020 | ESP | 24 | -0.0155 | [-0.0427, +0.0027] | -2.8 | 0.37 | [1, 132] | tie |
| 50 | paris_marathon 2022 | FRA | 56 | -0.0150 | [-0.0323, -0.0046] | -2.7 | 0.27 | [3, 95] | tie |
