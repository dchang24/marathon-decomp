# Top-N slowest / fastest race editions by v_j

Estimand: v_j under the beta=0 (bundling) gauge = race effect RELATIVE
TO CONTEMPORANEOUS RACES; it bundles course + that day's weather +
field conditions. Selection is by bootstrap rank stability P(top-N)
(headline P>=0.5, tie P>=0.25), which demotes small-field races whose
extreme point v_j is noise (winner's curse). min@3:00 =
180*(exp(v_j)-1) = minutes vs the average race for a 3:00:00 runner.

slice: ALL_M_14-25_mrc2; min_n = 0; n_top = 20; model = full_nu8p00

### ALL_M_14-25_mrc2 -- slowest 20 (headline P>=0.5, tie P>=0.25)

| rank | race | country | n_j | v_j | v 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|
| 1 | milton_keynes_marathon 2018 | GBR | 398 | +0.1089 | [+0.1004, +0.1177] | +20.7 | 1.00 | [1, 1] | headline |
| 2 | nyc_marathon 2022 | USA | 13714 | +0.0844 | [+0.0825, +0.0855] | +15.9 | 1.00 | [2, 2] | headline |
| 3 | stockholm_marathon 2024 | SWE | 5198 | +0.0767 | [+0.0736, +0.0795] | +14.3 | 1.00 | [3, 3] | headline |
| 4 | london_marathon 2018 | GBR | 5885 | +0.0636 | [+0.0611, +0.0654] | +11.8 | 1.00 | [4, 4] | headline |
| 5 | sydney_marathon 2023 | AUS | 1962 | +0.0481 | [+0.0434, +0.0516] | +8.9 | 1.00 | [5, 8] | headline |
| 6 | boston_marathon 2024 | USA | 2579 | +0.0457 | [+0.0434, +0.0479] | +8.4 | 1.00 | [5, 8] | headline |
| 7 | athens_marathon 2025 | GRE | 669 | +0.0435 | [+0.0376, +0.0492] | +8.0 | 1.00 | [5, 12] | headline |
| 8 | boston_marathon 2017 | USA | 1658 | +0.0429 | [+0.0404, +0.0462] | +7.9 | 1.00 | [6, 10] | headline |
| 9 | boston_marathon 2018 | USA | 1893 | +0.0376 | [+0.0341, +0.0405] | +6.9 | 1.00 | [9, 16] | headline |
| 10 | milton_keynes_marathon 2016 | GBR | 380 | +0.0375 | [+0.0291, +0.0468] | +6.9 | 0.97 | [6, 21] | headline |
| 11 | melbourne_marathon 2016 | AUS | 991 | +0.0369 | [+0.0231, +0.0445] | +6.8 | 0.96 | [7, 35] | headline |
| 12 | berlin_marathon 2025 | GER | 3128 | +0.0365 | [+0.0350, +0.0389] | +6.7 | 1.00 | [9, 15] | headline |
| 13 | stockholm_marathon 2022 | SWE | 3869 | +0.0354 | [+0.0325, +0.0383] | +6.5 | 1.00 | [10, 18] | headline |
| 14 | nyc_marathon 2023 | USA | 14374 | +0.0350 | [+0.0336, +0.0362] | +6.4 | 1.00 | [11, 18] | headline |
| 15 | chicago_marathon 2021 | USA | 4323 | +0.0333 | [+0.0305, +0.0355] | +6.1 | 1.00 | [13, 20] | headline |
| 16 | milton_keynes_marathon 2023 | GBR | 297 | +0.0325 | [+0.0254, +0.0391] | +6.0 | 0.87 | [10, 29] | headline |
| 17 | melbourne_marathon 2018 | AUS | 1184 | +0.0313 | [+0.0259, +0.0371] | +5.7 | 0.74 | [11, 27] | headline |
| 18 | boston_marathon 2016 | USA | 1273 | +0.0305 | [+0.0279, +0.0345] | +5.6 | 0.81 | [15, 23] | headline |
| 19 | vienna_marathon 2021 | AUT | 232 | +0.0274 | [+0.0224, +0.0344] | +5.0 | 0.31 | [14, 38] | tie |
| 24 | stockholm_marathon 2018 | SWE | 455 | +0.0264 | [+0.0216, +0.0318] | +4.8 | 0.25 | [17, 40] | tie |
| 28 | belfast_marathon 2017 | GBR | 74 | +0.0237 | [+0.0109, +0.0386] | +4.3 | 0.25 | [11, 90] | tie |
| 32 | lisbon_marathon 2023 | POR | 52 | +0.0228 | [+0.0110, +0.0362] | +4.1 | 0.25 | [12, 89] | tie |
| 34 | lisbon_marathon 2025 | POR | 32 | +0.0225 | [+0.0063, +0.0475] | +4.1 | 0.25 | [5, 114] | tie |

### ALL_M_14-25_mrc2 -- fastest 20 (headline P>=0.5, tie P>=0.25)

| rank | race | country | n_j | v_j | v 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|
| 1 | chicago_marathon 2014 | USA | 5269 | -0.0355 | [-0.0368, -0.0324] | -6.3 | 1.00 | [1, 4] | headline |
| 2 | eindhoven_marathon 2025 | NED | 795 | -0.0335 | [-0.0365, -0.0290] | -5.9 | 1.00 | [1, 7] | headline |
| 3 | amsterdam_marathon 2025 | NED | 4046 | -0.0333 | [-0.0353, -0.0304] | -5.9 | 1.00 | [1, 6] | headline |
| 4 | berlin_marathon 2014 | GER | 3196 | -0.0327 | [-0.0346, -0.0299] | -5.8 | 1.00 | [2, 6] | headline |
| 5 | seville_marathon 2019 | ESP | 394 | -0.0276 | [-0.0343, -0.0230] | -4.9 | 0.94 | [2, 24] | headline |
| 6 | eindhoven_marathon 2024 | NED | 1155 | -0.0270 | [-0.0299, -0.0236] | -4.8 | 0.98 | [5, 19] | headline |
| 7 | florence_marathon 2025 | ITA | 2967 | -0.0266 | [-0.0295, -0.0241] | -4.7 | 0.99 | [6, 17] | headline |
| 8 | hamburg_marathon 2019 | GER | 2873 | -0.0263 | [-0.0282, -0.0231] | -4.7 | 0.90 | [6, 25] | headline |
| 9 | seville_marathon 2018 | ESP | 310 | -0.0263 | [-0.0323, -0.0208] | -4.7 | 0.87 | [3, 30] | headline |
| 10 | stockholm_marathon 2015 | SWE | 1747 | -0.0262 | [-0.0294, -0.0205] | -4.6 | 0.73 | [6, 35] | headline |
| 11 | zurich_marathon 2025 | SUI | 982 | -0.0253 | [-0.0289, -0.0215] | -4.5 | 0.77 | [6, 29] | headline |
| 12 | berlin_marathon 2015 | GER | 4501 | -0.0253 | [-0.0263, -0.0225] | -4.5 | 0.85 | [8, 25] | headline |
| 13 | frankfurt_marathon 2015 | GER | 2318 | -0.0252 | [-0.0273, -0.0228] | -4.5 | 0.90 | [7, 23] | headline |
| 14 | prague_marathon 2014 | CZE | 239 | -0.0247 | [-0.0346, -0.0142] | -4.4 | 0.50 | [1, 87] | headline |
| 15 | amsterdam_marathon 2023 | NED | 3038 | -0.0243 | [-0.0262, -0.0221] | -4.3 | 0.67 | [9, 28] | headline |
| 16 | london_marathon 2015 | GBR | 5597 | -0.0242 | [-0.0261, -0.0214] | -4.3 | 0.53 | [13, 29] | headline |
| 17 | newport_marathon 2023 | GBR | 352 | -0.0240 | [-0.0297, -0.0174] | -4.3 | 0.56 | [6, 55] | headline |
| 18 | rotterdam_marathon 2023 | NED | 6091 | -0.0239 | [-0.0259, -0.0222] | -4.2 | 0.67 | [11, 28] | headline |
| 19 | valencia_marathon 2023 | ESP | 8820 | -0.0234 | [-0.0247, -0.0220] | -4.2 | 0.38 | [15, 29] | tie |
| 20 | paris_marathon 2023 | FRA | 10006 | -0.0232 | [-0.0267, -0.0209] | -4.1 | 0.51 | [10, 34] | headline |
| 21 | copenhagen_marathon 2014 | DEN | 1266 | -0.0232 | [-0.0276, -0.0189] | -4.1 | 0.40 | [7, 46] | tie |
| 22 | zurich_marathon 2019 | SUI | 553 | -0.0226 | [-0.0292, -0.0167] | -4.0 | 0.28 | [6, 61] | tie |
| 23 | dublin_marathon 2019 | IRL | 150 | -0.0225 | [-0.0293, -0.0143] | -4.0 | 0.30 | [5, 84] | tie |
| 24 | prague_marathon 2019 | CZE | 329 | -0.0223 | [-0.0295, -0.0159] | -4.0 | 0.39 | [5, 71] | tie |
| 25 | valencia_marathon 2020 | ESP | 51 | -0.0217 | [-0.0282, -0.0154] | -3.9 | 0.35 | [7, 75] | tie |
| 26 | tokyo_marathon 2022 | JPN | 176 | -0.0215 | [-0.0285, -0.0156] | -3.8 | 0.33 | [6, 71] | tie |
| 30 | bostonUK_marathon 2023 | GBR | 159 | -0.0208 | [-0.0325, -0.0140] | -3.7 | 0.39 | [3, 87] | tie |
| 32 | newport_marathon 2019 | GBR | 352 | -0.0206 | [-0.0330, -0.0142] | -3.7 | 0.39 | [4, 87] | tie |
| 90 | san_sebastian_marathon 2023 | ESP | 29 | -0.0138 | [-0.0470, +0.0183] | -2.5 | 0.31 | [1, 291] | tie |
