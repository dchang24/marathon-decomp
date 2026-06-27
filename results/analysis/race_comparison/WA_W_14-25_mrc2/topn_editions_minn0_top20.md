# Top-N slowest / fastest race editions by v_j

Estimand: v_j under the beta=0 (bundling) gauge = race effect RELATIVE
TO CONTEMPORANEOUS RACES; it bundles course + that day's weather +
field conditions. Selection is by bootstrap rank stability P(top-N)
(headline P>=0.5, tie P>=0.25), which demotes small-field races whose
extreme point v_j is noise (winner's curse). min@3:00 =
180*(exp(v_j)-1) = minutes vs the average race for a 3:00:00 runner.

slice: WA_W_14-25_mrc2; min_n = 0; n_top = 20; model = full_nu8p00

### WA_W_14-25_mrc2 -- slowest 20 (headline P>=0.5, tie P>=0.25)

| rank | race | country | n_j | v_j | v 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|
| 1 | stockholm_marathon 2024 | SWE | 31 | +0.0559 | [+0.0442, +0.0679] | +10.3 | 1.00 | [1, 3] | headline |
| 2 | nyc_marathon 2022 | USA | 142 | +0.0467 | [+0.0381, +0.0547] | +8.6 | 1.00 | [1, 4] | headline |
| 3 | london_marathon 2018 | GBR | 150 | +0.0467 | [+0.0367, +0.0573] | +8.6 | 1.00 | [1, 5] | headline |
| 4 | boston_marathon 2018 | USA | 117 | +0.0374 | [+0.0313, +0.0443] | +6.9 | 1.00 | [3, 7] | headline |
| 5 | nyc_marathon 2025 | USA | 115 | +0.0331 | [+0.0230, +0.0427] | +6.1 | 1.00 | [3, 13] | headline |
| 6 | nyc_marathon 2024 | USA | 133 | +0.0279 | [+0.0204, +0.0342] | +5.1 | 0.98 | [5, 19] | headline |
| 7 | stockholm_marathon 2022 | SWE | 26 | +0.0271 | [+0.0078, +0.0441] | +4.9 | 0.85 | [3, 38] | headline |
| 8 | berlin_marathon 2025 | GER | 144 | +0.0270 | [+0.0201, +0.0350] | +4.9 | 1.00 | [5, 18] | headline |
| 9 | boston_marathon 2017 | USA | 147 | +0.0262 | [+0.0202, +0.0312] | +4.8 | 1.00 | [6, 18] | headline |
| 10 | nyc_marathon 2023 | USA | 110 | +0.0257 | [+0.0179, +0.0358] | +4.7 | 0.99 | [4, 20] | headline |
| 11 | nyc_marathon 2017 | USA | 131 | +0.0255 | [+0.0185, +0.0316] | +4.6 | 0.98 | [7, 18] | headline |
| 12 | boston_marathon 2016 | USA | 122 | +0.0238 | [+0.0173, +0.0321] | +4.3 | 0.96 | [6, 21] | headline |
| 13 | nyc_marathon 2014 | USA | 83 | +0.0204 | [+0.0105, +0.0312] | +3.7 | 0.77 | [6, 30] | headline |
| 14 | stockholm_marathon 2023 | SWE | 26 | +0.0204 | [+0.0074, +0.0334] | +3.7 | 0.66 | [5, 39] | headline |
| 15 | nyc_marathon 2015 | USA | 106 | +0.0202 | [+0.0125, +0.0276] | +3.7 | 0.73 | [8, 28] | headline |
| 16 | manchester_marathon 2021 | GBR | 30 | +0.0200 | [+0.0043, +0.0416] | +3.6 | 0.57 | [4, 48] | headline |
| 17 | stockholm_marathon 2025 | SWE | 32 | +0.0198 | [+0.0109, +0.0290] | +3.6 | 0.73 | [8, 31] | headline |
| 18 | nyc_marathon 2016 | USA | 115 | +0.0193 | [+0.0115, +0.0264] | +3.5 | 0.65 | [10, 33] | headline |
| 19 | sydney_marathon 2024 | AUS | 45 | +0.0185 | [+0.0064, +0.0307] | +3.4 | 0.63 | [6, 44] | headline |
| 20 | boston_marathon 2024 | USA | 136 | +0.0183 | [+0.0118, +0.0236] | +3.3 | 0.53 | [12, 31] | headline |
| 21 | milan_marathon 2025 | ITA | 22 | +0.0181 | [-0.0001, +0.0419] | +3.3 | 0.50 | [4, 65] | headline |
| 22 | boston_marathon 2019 | USA | 209 | +0.0176 | [+0.0110, +0.0226] | +3.2 | 0.49 | [14, 30] | tie |
| 23 | rotterdam_marathon 2018 | NED | 36 | +0.0148 | [+0.0025, +0.0254] | +2.7 | 0.27 | [10, 57] | tie |

### WA_W_14-25_mrc2 -- fastest 20 (headline P>=0.5, tie P>=0.25)

| rank | race | country | n_j | v_j | v 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|
| 1 | dublin_marathon 2025 | IRL | 24 | -0.0259 | [-0.0355, -0.0128] | -4.6 | 0.93 | [1, 30] | headline |
| 2 | paris_marathon 2015 | FRA | 20 | -0.0252 | [-0.0345, -0.0113] | -4.5 | 0.95 | [1, 36] | headline |
| 3 | amsterdam_marathon 2025 | NED | 50 | -0.0223 | [-0.0336, -0.0125] | -4.0 | 0.94 | [1, 28] | headline |
| 4 | dublin_marathon 2024 | IRL | 33 | -0.0219 | [-0.0324, -0.0094] | -3.9 | 0.91 | [1, 45] | headline |
| 5 | tokyo_marathon 2015 | JPN | 53 | -0.0219 | [-0.0332, -0.0119] | -3.9 | 0.95 | [1, 31] | headline |
| 6 | berlin_marathon 2014 | GER | 43 | -0.0202 | [-0.0300, -0.0095] | -3.6 | 0.77 | [1, 43] | headline |
| 7 | seville_marathon 2020 | ESP | 32 | -0.0201 | [-0.0336, -0.0062] | -3.6 | 0.81 | [1, 58] | headline |
| 8 | dublin_marathon 2019 | IRL | 39 | -0.0196 | [-0.0300, -0.0044] | -3.5 | 0.75 | [1, 67] | headline |
| 9 | valencia_marathon 2023 | ESP | 267 | -0.0194 | [-0.0228, -0.0159] | -3.5 | 0.99 | [4, 19] | headline |
| 10 | valencia_marathon 2022 | ESP | 157 | -0.0184 | [-0.0223, -0.0136] | -3.3 | 0.85 | [5, 25] | headline |
| 11 | manchester_marathon 2022 | GBR | 78 | -0.0174 | [-0.0251, -0.0102] | -3.1 | 0.70 | [4, 44] | headline |
| 12 | tokyo_marathon 2016 | JPN | 61 | -0.0169 | [-0.0259, -0.0092] | -3.0 | 0.54 | [4, 47] | headline |
| 13 | seville_marathon 2022 | ESP | 51 | -0.0164 | [-0.0231, -0.0076] | -2.9 | 0.52 | [5, 52] | headline |
| 14 | barcelona_marathon 2024 | ESP | 27 | -0.0160 | [-0.0248, -0.0046] | -2.9 | 0.51 | [4, 68] | headline |
| 15 | rotterdam_marathon 2023 | NED | 74 | -0.0160 | [-0.0213, -0.0100] | -2.9 | 0.65 | [6, 39] | headline |
| 16 | chicago_marathon 2014 | USA | 111 | -0.0159 | [-0.0218, -0.0082] | -2.8 | 0.52 | [6, 50] | headline |
| 17 | frankfurt_marathon 2016 | GER | 51 | -0.0158 | [-0.0249, -0.0090] | -2.8 | 0.55 | [2, 48] | headline |
| 18 | amsterdam_marathon 2023 | NED | 53 | -0.0149 | [-0.0211, -0.0078] | -2.7 | 0.42 | [7, 51] | tie |
| 19 | barcelona_marathon 2025 | ESP | 52 | -0.0145 | [-0.0233, -0.0069] | -2.6 | 0.45 | [4, 58] | tie |
| 20 | malaga_marathon 2025 | ESP | 27 | -0.0141 | [-0.0262, -0.0036] | -2.5 | 0.42 | [4, 70] | tie |
| 21 | berlin_marathon 2022 | GER | 144 | -0.0140 | [-0.0188, -0.0083] | -2.5 | 0.30 | [11, 50] | tie |
| 22 | chicago_marathon 2023 | USA | 235 | -0.0139 | [-0.0185, -0.0093] | -2.5 | 0.32 | [12, 46] | tie |
| 23 | seville_marathon 2019 | ESP | 22 | -0.0135 | [-0.0259, +0.0020] | -2.4 | 0.39 | [2, 93] | tie |
| 24 | chicago_marathon 2019 | USA | 195 | -0.0130 | [-0.0198, -0.0070] | -2.3 | 0.33 | [9, 58] | tie |
| 25 | hamburg_marathon 2025 | GER | 23 | -0.0129 | [-0.0254, +0.0003] | -2.3 | 0.28 | [3, 86] | tie |
| 26 | seville_marathon 2024 | ESP | 85 | -0.0128 | [-0.0196, -0.0064] | -2.3 | 0.29 | [9, 56] | tie |
| 27 | manchester_marathon 2023 | GBR | 34 | -0.0128 | [-0.0269, +0.0007] | -2.3 | 0.44 | [2, 89] | tie |
| 28 | frankfurt_marathon 2015 | GER | 37 | -0.0126 | [-0.0220, -0.0042] | -2.3 | 0.33 | [6, 68] | tie |
| 44 | paris_marathon 2019 | FRA | 24 | -0.0100 | [-0.0233, +0.0018] | -1.8 | 0.26 | [4, 93] | tie |
