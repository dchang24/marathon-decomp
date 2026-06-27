# Top-N slowest / fastest race editions by v_j

Estimand: v_j under the beta=0 (bundling) gauge = race effect RELATIVE
TO CONTEMPORANEOUS RACES; it bundles course + that day's weather +
field conditions. Selection is by bootstrap rank stability P(top-N)
(headline P>=0.5, tie P>=0.25), which demotes small-field races whose
extreme point v_j is noise (winner's curse). min@3:00 =
180*(exp(v_j)-1) = minutes vs the average race for a 3:00:00 runner.

slice: ALL_W_14-25_mrc2; min_n = 0; n_top = 20; model = full_nu8p00

### ALL_W_14-25_mrc2 -- slowest 20 (headline P>=0.5, tie P>=0.25)

| rank | race | country | n_j | v_j | v 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|
| 1 | milton_keynes_marathon 2018 | GBR | 137 | +0.1141 | [+0.1012, +0.1264] | +21.8 | 1.00 | [1, 1] | headline |
| 2 | nyc_marathon 2022 | USA | 9979 | +0.0684 | [+0.0665, +0.0704] | +12.7 | 1.00 | [2, 3] | headline |
| 3 | stockholm_marathon 2024 | SWE | 1786 | +0.0658 | [+0.0624, +0.0701] | +12.2 | 1.00 | [2, 3] | headline |
| 4 | london_marathon 2018 | GBR | 3277 | +0.0608 | [+0.0579, +0.0632] | +11.3 | 1.00 | [4, 5] | headline |
| 5 | stockholm_marathon 2018 | SWE | 1589 | +0.0560 | [+0.0525, +0.0617] | +10.4 | 1.00 | [4, 6] | headline |
| 6 | sydney_marathon 2023 | AUS | 843 | +0.0522 | [+0.0479, +0.0572] | +9.6 | 1.00 | [5, 7] | headline |
| 7 | athens_marathon 2025 | GRE | 136 | +0.0411 | [+0.0296, +0.0502] | +7.6 | 0.99 | [6, 14] | headline |
| 8 | berlin_marathon 2025 | GER | 1537 | +0.0408 | [+0.0388, +0.0442] | +7.5 | 1.00 | [7, 9] | headline |
| 9 | boston_marathon 2024 | USA | 1476 | +0.0360 | [+0.0332, +0.0389] | +6.6 | 1.00 | [8, 13] | headline |
| 10 | nyc_marathon 2023 | USA | 10346 | +0.0344 | [+0.0324, +0.0357] | +6.3 | 1.00 | [9, 14] | headline |
| 11 | vienna_marathon 2018 | AUT | 106 | +0.0331 | [+0.0244, +0.0410] | +6.1 | 0.90 | [8, 25] | headline |
| 12 | madrid_marathon 2023 | ESP | 43 | +0.0310 | [+0.0152, +0.0433] | +5.7 | 0.75 | [7, 60] | headline |
| 13 | belfast_marathon 2017 | GBR | 27 | +0.0295 | [+0.0066, +0.0464] | +5.4 | 0.60 | [7, 102] | headline |
| 14 | boston_marathon 2018 | USA | 1114 | +0.0279 | [+0.0247, +0.0310] | +5.1 | 0.77 | [13, 25] | headline |
| 15 | venice_marathon 2023 | ITA | 329 | +0.0278 | [+0.0222, +0.0330] | +5.1 | 0.56 | [12, 33] | headline |
| 16 | chicago_marathon 2021 | USA | 3696 | +0.0277 | [+0.0257, +0.0302] | +5.1 | 0.82 | [14, 23] | headline |
| 17 | rome_marathon 2021 | ITA | 396 | +0.0273 | [+0.0202, +0.0332] | +5.0 | 0.70 | [11, 39] | headline |
| 18 | athens_marathon 2024 | GRE | 163 | +0.0266 | [+0.0184, +0.0338] | +4.8 | 0.44 | [12, 46] | tie |
| 19 | eindhoven_marathon 2018 | NED | 106 | +0.0264 | [+0.0169, +0.0337] | +4.8 | 0.52 | [12, 56] | headline |
| 20 | milton_keynes_marathon 2023 | GBR | 143 | +0.0256 | [+0.0149, +0.0356] | +4.7 | 0.38 | [11, 69] | tie |
| 21 | edinburgh_marathon 2017 | GBR | 198 | +0.0247 | [+0.0151, +0.0348] | +4.5 | 0.38 | [11, 59] | tie |
| 23 | prague_marathon 2022 | CZE | 74 | +0.0243 | [+0.0100, +0.0350] | +4.4 | 0.31 | [10, 81] | tie |
| 24 | melbourne_marathon 2018 | AUS | 441 | +0.0236 | [+0.0151, +0.0356] | +4.3 | 0.38 | [10, 63] | tie |
| 26 | vienna_marathon 2021 | AUT | 43 | +0.0235 | [+0.0121, +0.0361] | +4.3 | 0.28 | [10, 83] | tie |
| 31 | belfast_marathon 2018 | GBR | 47 | +0.0227 | [+0.0026, +0.0358] | +4.1 | 0.28 | [11, 130] | tie |
| 51 | madrid_marathon 2021 | ESP | 25 | +0.0172 | [-0.0041, +0.0391] | +3.1 | 0.25 | [9, 181] | tie |

### ALL_W_14-25_mrc2 -- fastest 20 (headline P>=0.5, tie P>=0.25)

| rank | race | country | n_j | v_j | v 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|
| 1 | eindhoven_marathon 2025 | NED | 158 | -0.0329 | [-0.0410, -0.0253] | -5.8 | 1.00 | [1, 13] | headline |
| 2 | chicago_marathon 2014 | USA | 3876 | -0.0311 | [-0.0333, -0.0268] | -5.5 | 1.00 | [1, 10] | headline |
| 3 | amsterdam_marathon 2025 | NED | 1086 | -0.0305 | [-0.0329, -0.0265] | -5.4 | 1.00 | [1, 10] | headline |
| 4 | belfast_marathon 2014 | GBR | 21 | -0.0264 | [-0.0566, +0.0024] | -4.7 | 0.62 | [1, 208] | headline |
| 5 | zurich_marathon 2025 | SUI | 204 | -0.0258 | [-0.0344, -0.0150] | -4.6 | 0.86 | [1, 62] | headline |
| 6 | eindhoven_marathon 2024 | NED | 234 | -0.0245 | [-0.0311, -0.0183] | -4.3 | 0.74 | [4, 38] | headline |
| 7 | berlin_marathon 2014 | GER | 892 | -0.0244 | [-0.0284, -0.0199] | -4.3 | 0.87 | [5, 33] | headline |
| 8 | dublin_marathon 2019 | IRL | 83 | -0.0244 | [-0.0315, -0.0168] | -4.3 | 0.66 | [2, 49] | headline |
| 9 | chester_marathon 2016 | GBR | 168 | -0.0238 | [-0.0350, -0.0144] | -4.2 | 0.67 | [1, 67] | headline |
| 10 | newport_marathon 2025 | GBR | 150 | -0.0238 | [-0.0317, -0.0156] | -4.2 | 0.60 | [2, 55] | headline |
| 11 | seville_marathon 2018 | ESP | 60 | -0.0237 | [-0.0391, -0.0088] | -4.2 | 0.59 | [1, 116] | headline |
| 12 | hannover_marathon 2025 | GER | 127 | -0.0236 | [-0.0306, -0.0179] | -4.2 | 0.73 | [3, 44] | headline |
| 13 | dublin_marathon 2025 | IRL | 56 | -0.0228 | [-0.0324, -0.0083] | -4.1 | 0.45 | [2, 125] | tie |
| 14 | chester_marathon 2025 | GBR | 208 | -0.0225 | [-0.0317, -0.0132] | -4.0 | 0.49 | [3, 81] | tie |
| 15 | melbourne_marathon 2014 | AUS | 306 | -0.0225 | [-0.0643, +0.0071] | -4.0 | 0.50 | [1, 238] | headline |
| 16 | florence_marathon 2024 | ITA | 741 | -0.0223 | [-0.0273, -0.0191] | -4.0 | 0.72 | [8, 36] | headline |
| 17 | rotterdam_marathon 2023 | NED | 1662 | -0.0217 | [-0.0247, -0.0187] | -3.9 | 0.33 | [12, 39] | tie |
| 18 | boston_marathon 2015 | USA | 874 | -0.0217 | [-0.0246, -0.0167] | -3.9 | 0.31 | [12, 50] | tie |
| 19 | copenhagen_marathon 2025 | DEN | 803 | -0.0216 | [-0.0246, -0.0171] | -3.8 | 0.31 | [13, 48] | tie |
| 20 | seville_marathon 2019 | ESP | 84 | -0.0214 | [-0.0309, -0.0132] | -3.8 | 0.55 | [3, 80] | headline |
| 21 | paris_marathon 2023 | FRA | 2579 | -0.0214 | [-0.0251, -0.0168] | -3.8 | 0.38 | [10, 49] | tie |
| 22 | malaga_marathon 2018 | ESP | 54 | -0.0211 | [-0.0302, -0.0074] | -3.8 | 0.43 | [3, 131] | tie |
| 23 | london_marathon 2015 | GBR | 2649 | -0.0210 | [-0.0236, -0.0184] | -3.7 | 0.29 | [14, 42] | tie |
| 25 | copenhagen_marathon 2014 | DEN | 313 | -0.0204 | [-0.0296, -0.0130] | -3.6 | 0.41 | [5, 79] | tie |
| 27 | dublin_marathon 2024 | IRL | 63 | -0.0195 | [-0.0290, -0.0048] | -3.5 | 0.32 | [6, 157] | tie |
| 28 | cape_town_marathon 2023 | RSA | 29 | -0.0194 | [-0.0332, -0.0097] | -3.5 | 0.39 | [2, 110] | tie |
| 30 | melbourne_marathon 2015 | AUS | 396 | -0.0193 | [-0.0616, +0.0102] | -3.4 | 0.30 | [1, 255] | tie |
| 40 | newport_marathon 2019 | GBR | 149 | -0.0176 | [-0.0292, -0.0080] | -3.1 | 0.25 | [6, 123] | tie |
| 48 | malaga_marathon 2021 | ESP | 36 | -0.0165 | [-0.0296, -0.0049] | -2.9 | 0.28 | [5, 158] | tie |
