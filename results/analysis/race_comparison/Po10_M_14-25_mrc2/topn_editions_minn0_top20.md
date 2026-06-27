# Top-N slowest / fastest race editions by v_j

Estimand: v_j under the beta=0 (bundling) gauge = race effect RELATIVE
TO CONTEMPORANEOUS RACES; it bundles course + that day's weather +
field conditions. Selection is by bootstrap rank stability P(top-N)
(headline P>=0.5, tie P>=0.25), which demotes small-field races whose
extreme point v_j is noise (winner's curse). min@3:00 =
180*(exp(v_j)-1) = minutes vs the average race for a 3:00:00 runner.

slice: Po10_M_14-25_mrc2; min_n = 0; n_top = 20; model = full_nu8p00

### Po10_M_14-25_mrc2 -- slowest 20 (headline P>=0.5, tie P>=0.25)

| rank | race | country | n_j | v_j | v 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|
| 1 | milton_keynes_marathon 2018 | GBR | 396 | +0.1091 | [+0.1000, +0.1168] | +20.7 | 1.00 | [1, 2] | headline |
| 2 | stockholm_marathon 2024 | SWE | 24 | +0.0915 | [+0.0593, +0.1136] | +17.2 | 1.00 | [1, 8] | headline |
| 3 | nyc_marathon 2022 | USA | 368 | +0.0911 | [+0.0846, +0.0973] | +17.2 | 1.00 | [2, 4] | headline |
| 4 | london_marathon 2018 | GBR | 4417 | +0.0672 | [+0.0647, +0.0690] | +12.5 | 1.00 | [4, 8] | headline |
| 5 | vienna_marathon 2018 | AUT | 24 | +0.0642 | [+0.0403, +0.0944] | +11.9 | 0.99 | [2, 16] | headline |
| 6 | athens_marathon 2025 | GRE | 40 | +0.0640 | [+0.0438, +0.0872] | +11.9 | 1.00 | [3, 14] | headline |
| 7 | boston_marathon 2017 | USA | 190 | +0.0585 | [+0.0496, +0.0677] | +10.8 | 1.00 | [5, 11] | headline |
| 8 | chicago_marathon 2021 | USA | 23 | +0.0548 | [+0.0095, +0.0849] | +10.1 | 0.85 | [3, 68] | headline |
| 9 | venice_marathon 2024 | ITA | 21 | +0.0522 | [+0.0274, +0.0855] | +9.6 | 0.87 | [3, 30] | headline |
| 10 | boston_marathon 2024 | USA | 403 | +0.0504 | [+0.0450, +0.0570] | +9.3 | 1.00 | [7, 13] | headline |
| 11 | stockholm_marathon 2022 | SWE | 30 | +0.0485 | [+0.0242, +0.0734] | +9.0 | 0.95 | [4, 31] | headline |
| 12 | vienna_marathon 2023 | AUT | 38 | +0.0477 | [+0.0249, +0.0765] | +8.8 | 0.86 | [4, 31] | headline |
| 13 | berlin_marathon 2025 | GER | 654 | +0.0457 | [+0.0421, +0.0511] | +8.4 | 1.00 | [8, 15] | headline |
| 14 | copenhagen_marathon 2018 | DEN | 54 | +0.0405 | [+0.0211, +0.0594] | +7.4 | 0.73 | [7, 35] | headline |
| 15 | copenhagen_marathon 2016 | DEN | 48 | +0.0400 | [+0.0164, +0.0572] | +7.3 | 0.72 | [7, 51] | headline |
| 16 | boston_marathon 2016 | USA | 131 | +0.0394 | [+0.0304, +0.0512] | +7.2 | 0.88 | [9, 24] | headline |
| 17 | boston_marathon 2018 | USA | 209 | +0.0379 | [+0.0301, +0.0470] | +6.9 | 0.84 | [11, 26] | headline |
| 18 | stockholm_marathon 2023 | SWE | 37 | +0.0355 | [+0.0216, +0.0508] | +6.5 | 0.56 | [11, 35] | headline |
| 19 | milton_keynes_marathon 2023 | GBR | 297 | +0.0335 | [+0.0259, +0.0409] | +6.1 | 0.43 | [14, 31] | tie |
| 20 | boston_marathon 2021 | USA | 29 | +0.0329 | [+0.0109, +0.0530] | +6.0 | 0.46 | [10, 65] | tie |
| 21 | milton_keynes_marathon 2016 | GBR | 377 | +0.0326 | [+0.0255, +0.0405] | +6.0 | 0.42 | [14, 29] | tie |
| 23 | munich_marathon 2023 | GER | 27 | +0.0303 | [+0.0048, +0.0514] | +5.5 | 0.33 | [10, 82] | tie |
| 29 | prague_marathon 2024 | CZE | 21 | +0.0248 | [+0.0013, +0.0459] | +4.5 | 0.29 | [12, 102] | tie |

### Po10_M_14-25_mrc2 -- fastest 20 (headline P>=0.5, tie P>=0.25)

| rank | race | country | n_j | v_j | v 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|
| 1 | hamburg_marathon 2019 | GER | 23 | -0.0435 | [-0.0593, -0.0309] | -7.7 | 0.98 | [1, 10] | headline |
| 2 | malaga_marathon 2019 | ESP | 41 | -0.0417 | [-0.0521, -0.0296] | -7.3 | 1.00 | [1, 13] | headline |
| 3 | munich_marathon 2025 | GER | 23 | -0.0383 | [-0.0803, -0.0013] | -6.8 | 0.70 | [1, 148] | headline |
| 4 | amsterdam_marathon 2025 | NED | 149 | -0.0368 | [-0.0427, -0.0287] | -6.5 | 0.99 | [1, 15] | headline |
| 5 | amsterdam_marathon 2023 | NED | 224 | -0.0333 | [-0.0393, -0.0281] | -5.9 | 1.00 | [2, 14] | headline |
| 6 | berlin_marathon 2014 | GER | 371 | -0.0317 | [-0.0369, -0.0233] | -5.6 | 0.95 | [4, 28] | headline |
| 7 | barcelona_marathon 2021 | ESP | 39 | -0.0314 | [-0.0465, -0.0205] | -5.6 | 0.85 | [2, 36] | headline |
| 8 | london_marathon 2015 | GBR | 4875 | -0.0303 | [-0.0318, -0.0259] | -5.4 | 1.00 | [7, 18] | headline |
| 9 | seville_marathon 2019 | ESP | 67 | -0.0287 | [-0.0389, -0.0175] | -5.1 | 0.87 | [3, 55] | headline |
| 10 | chicago_marathon 2014 | USA | 90 | -0.0280 | [-0.0432, -0.0144] | -5.0 | 0.61 | [3, 71] | headline |
| 11 | copenhagen_marathon 2025 | DEN | 84 | -0.0278 | [-0.0350, -0.0167] | -4.9 | 0.50 | [6, 57] | headline |
| 12 | valencia_marathon 2023 | ESP | 587 | -0.0266 | [-0.0307, -0.0230] | -4.7 | 0.75 | [9, 28] | headline |
| 13 | dublin_marathon 2019 | IRL | 117 | -0.0264 | [-0.0359, -0.0173] | -4.7 | 0.62 | [4, 57] | headline |
| 14 | rotterdam_marathon 2023 | NED | 35 | -0.0255 | [-0.0380, -0.0143] | -4.5 | 0.57 | [3, 76] | headline |
| 15 | paris_marathon 2019 | FRA | 172 | -0.0254 | [-0.0348, -0.0161] | -4.5 | 0.47 | [6, 66] | tie |
| 16 | frankfurt_marathon 2015 | GER | 111 | -0.0252 | [-0.0362, -0.0123] | -4.5 | 0.45 | [5, 81] | tie |
| 17 | tokyo_marathon 2015 | JPN | 38 | -0.0251 | [-0.0408, -0.0107] | -4.5 | 0.34 | [3, 93] | tie |
| 18 | boston_marathon 2015 | USA | 116 | -0.0250 | [-0.0304, -0.0164] | -4.4 | 0.34 | [9, 59] | tie |
| 19 | barcelona_marathon 2025 | ESP | 221 | -0.0247 | [-0.0313, -0.0161] | -4.4 | 0.49 | [9, 62] | tie |
| 20 | paris_marathon 2025 | FRA | 66 | -0.0246 | [-0.0354, -0.0115] | -4.4 | 0.38 | [4, 93] | tie |
| 21 | vienna_marathon 2025 | AUT | 58 | -0.0245 | [-0.0403, -0.0094] | -4.4 | 0.42 | [3, 104] | tie |
| 22 | newport_marathon 2023 | GBR | 350 | -0.0235 | [-0.0306, -0.0174] | -4.2 | 0.39 | [9, 56] | tie |
| 23 | san_sebastian_marathon 2025 | ESP | 34 | -0.0235 | [-0.0422, -0.0022] | -4.2 | 0.47 | [3, 135] | tie |
| 24 | seville_marathon 2020 | ESP | 90 | -0.0235 | [-0.0332, -0.0129] | -4.2 | 0.34 | [6, 84] | tie |
| 26 | paris_marathon 2023 | FRA | 161 | -0.0228 | [-0.0333, -0.0150] | -4.1 | 0.42 | [7, 69] | tie |
| 31 | newport_marathon 2019 | GBR | 351 | -0.0211 | [-0.0315, -0.0129] | -3.8 | 0.25 | [9, 80] | tie |
