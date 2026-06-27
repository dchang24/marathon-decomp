# Top-N slowest / fastest race editions by v_j

Estimand: v_j under the beta=0 (bundling) gauge = race effect RELATIVE
TO CONTEMPORANEOUS RACES; it bundles course + that day's weather +
field conditions. Selection is by bootstrap rank stability P(top-N)
(headline P>=0.5, tie P>=0.25), which demotes small-field races whose
extreme point v_j is noise (winner's curse). min@3:00 =
180*(exp(v_j)-1) = minutes vs the average race for a 3:00:00 runner.

slice: ALL_B_14-25_mrc2; min_n = 0; n_top = 20; model = full_nu8p00

### ALL_B_14-25_mrc2 -- slowest 20 (headline P>=0.5, tie P>=0.25)

| rank | race | country | n_j | v_j | v 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|
| 1 | milton_keynes_marathon 2018 | GBR | 535 | +0.1097 | [+0.1030, +0.1173] | +20.9 | 1.00 | [1, 1] | headline |
| 2 | nyc_marathon 2022 | USA | 23693 | +0.0772 | [+0.0758, +0.0781] | +14.4 | 1.00 | [2, 2] | headline |
| 3 | stockholm_marathon 2024 | SWE | 6984 | +0.0735 | [+0.0705, +0.0761] | +13.7 | 1.00 | [3, 3] | headline |
| 4 | london_marathon 2018 | GBR | 9162 | +0.0620 | [+0.0598, +0.0631] | +11.5 | 1.00 | [4, 4] | headline |
| 5 | sydney_marathon 2023 | AUS | 2805 | +0.0489 | [+0.0461, +0.0519] | +9.0 | 1.00 | [5, 7] | headline |
| 6 | stockholm_marathon 2018 | SWE | 2044 | +0.0476 | [+0.0446, +0.0506] | +8.8 | 1.00 | [5, 8] | headline |
| 7 | athens_marathon 2025 | GRE | 805 | +0.0432 | [+0.0386, +0.0493] | +7.9 | 1.00 | [5, 10] | headline |
| 8 | boston_marathon 2024 | USA | 4055 | +0.0417 | [+0.0400, +0.0439] | +7.7 | 1.00 | [7, 10] | headline |
| 9 | HK_marathon 2018 | HKG | 21 | +0.0412 | [+0.0321, +0.0508] | +7.6 | 1.00 | [5, 14] | headline |
| 10 | berlin_marathon 2025 | GER | 4665 | +0.0380 | [+0.0370, +0.0401] | +7.0 | 1.00 | [9, 11] | headline |
| 11 | boston_marathon 2017 | USA | 2731 | +0.0348 | [+0.0326, +0.0371] | +6.4 | 1.00 | [10, 16] | headline |
| 12 | nyc_marathon 2023 | USA | 24720 | +0.0346 | [+0.0333, +0.0356] | +6.3 | 1.00 | [11, 15] | headline |
| 13 | boston_marathon 2018 | USA | 3007 | +0.0334 | [+0.0315, +0.0362] | +6.1 | 1.00 | [11, 17] | headline |
| 14 | milton_keynes_marathon 2016 | GBR | 535 | +0.0325 | [+0.0251, +0.0393] | +5.9 | 0.90 | [10, 26] | headline |
| 15 | stockholm_marathon 2022 | SWE | 5243 | +0.0319 | [+0.0292, +0.0336] | +5.8 | 1.00 | [13, 19] | headline |
| 16 | chicago_marathon 2021 | USA | 8019 | +0.0304 | [+0.0287, +0.0319] | +5.6 | 0.98 | [15, 20] | headline |
| 17 | milton_keynes_marathon 2023 | GBR | 440 | +0.0303 | [+0.0233, +0.0358] | +5.5 | 0.74 | [11, 31] | headline |
| 18 | melbourne_marathon 2018 | AUS | 1625 | +0.0291 | [+0.0234, +0.0333] | +5.3 | 0.68 | [14, 31] | headline |
| 19 | venice_marathon 2023 | ITA | 1798 | +0.0273 | [+0.0242, +0.0308] | +5.0 | 0.40 | [16, 28] | tie |
| 20 | vienna_marathon 2021 | AUT | 275 | +0.0272 | [+0.0220, +0.0330] | +5.0 | 0.46 | [15, 35] | tie |
| 23 | lisbon_marathon 2025 | POR | 41 | +0.0251 | [+0.0105, +0.0448] | +4.6 | 0.48 | [7, 89] | tie |
| 26 | belfast_marathon 2017 | GBR | 101 | +0.0245 | [+0.0143, +0.0349] | +4.5 | 0.30 | [12, 69] | tie |

### ALL_B_14-25_mrc2 -- fastest 20 (headline P>=0.5, tie P>=0.25)

| rank | race | country | n_j | v_j | v 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|
| 1 | chicago_marathon 2014 | USA | 9145 | -0.0348 | [-0.0357, -0.0319] | -6.2 | 1.00 | [1, 4] | headline |
| 2 | eindhoven_marathon 2025 | NED | 953 | -0.0340 | [-0.0377, -0.0300] | -6.0 | 1.00 | [1, 5] | headline |
| 3 | amsterdam_marathon 2025 | NED | 5132 | -0.0330 | [-0.0342, -0.0306] | -5.8 | 1.00 | [1, 5] | headline |
| 4 | berlin_marathon 2014 | GER | 4088 | -0.0307 | [-0.0323, -0.0282] | -5.4 | 1.00 | [3, 7] | headline |
| 5 | eindhoven_marathon 2024 | NED | 1389 | -0.0270 | [-0.0297, -0.0240] | -4.8 | 1.00 | [5, 13] | headline |
| 6 | seville_marathon 2019 | ESP | 478 | -0.0263 | [-0.0315, -0.0216] | -4.7 | 0.97 | [3, 21] | headline |
| 7 | zurich_marathon 2025 | SUI | 1186 | -0.0255 | [-0.0292, -0.0218] | -4.5 | 0.91 | [5, 24] | headline |
| 8 | florence_marathon 2025 | ITA | 3669 | -0.0253 | [-0.0272, -0.0225] | -4.5 | 0.98 | [6, 20] | headline |
| 9 | seville_marathon 2018 | ESP | 372 | -0.0250 | [-0.0311, -0.0181] | -4.5 | 0.81 | [4, 40] | headline |
| 10 | hamburg_marathon 2019 | GER | 3716 | -0.0245 | [-0.0263, -0.0218] | -4.4 | 0.92 | [6, 24] | headline |
| 11 | berlin_marathon 2015 | GER | 5843 | -0.0238 | [-0.0248, -0.0218] | -4.2 | 0.90 | [10, 26] | headline |
| 12 | london_marathon 2015 | GBR | 8246 | -0.0238 | [-0.0252, -0.0221] | -4.2 | 0.95 | [9, 23] | headline |
| 13 | dublin_marathon 2019 | IRL | 233 | -0.0236 | [-0.0303, -0.0169] | -4.2 | 0.60 | [5, 51] | headline |
| 14 | rotterdam_marathon 2023 | NED | 7753 | -0.0236 | [-0.0253, -0.0223] | -4.2 | 0.95 | [9, 21] | headline |
| 15 | paris_marathon 2023 | FRA | 12585 | -0.0231 | [-0.0253, -0.0213] | -4.1 | 0.74 | [9, 26] | headline |
| 16 | amsterdam_marathon 2023 | NED | 3896 | -0.0231 | [-0.0244, -0.0210] | -4.1 | 0.79 | [12, 27] | headline |
| 17 | frankfurt_marathon 2015 | GER | 2923 | -0.0230 | [-0.0248, -0.0209] | -4.1 | 0.63 | [10, 28] | headline |
| 18 | san_sebastian_marathon 2023 | ESP | 33 | -0.0227 | [-0.0482, +0.0033] | -4.0 | 0.54 | [1, 213] | headline |
| 19 | valencia_marathon 2023 | ESP | 10350 | -0.0224 | [-0.0236, -0.0211] | -4.0 | 0.43 | [14, 27] | tie |
| 20 | stockholm_marathon 2015 | SWE | 2280 | -0.0222 | [-0.0253, -0.0179] | -3.9 | 0.32 | [8, 44] | tie |
| 21 | newport_marathon 2023 | GBR | 481 | -0.0222 | [-0.0285, -0.0172] | -3.9 | 0.55 | [5, 46] | headline |
| 22 | copenhagen_marathon 2014 | DEN | 1579 | -0.0220 | [-0.0261, -0.0189] | -3.9 | 0.41 | [7, 37] | tie |
| 23 | tokyo_marathon 2015 | JPN | 1182 | -0.0217 | [-0.0250, -0.0169] | -3.9 | 0.26 | [10, 54] | tie |
| 29 | newport_marathon 2019 | GBR | 501 | -0.0200 | [-0.0268, -0.0139] | -3.6 | 0.27 | [8, 81] | tie |
| 30 | prague_marathon 2019 | CZE | 404 | -0.0198 | [-0.0253, -0.0144] | -3.5 | 0.32 | [9, 78] | tie |
