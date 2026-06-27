# Top-N slowest / fastest race editions by v_j

Estimand: v_j under the beta=0 (bundling) gauge = race effect RELATIVE
TO CONTEMPORANEOUS RACES; it bundles course + that day's weather +
field conditions. Selection is by bootstrap rank stability P(top-N)
(headline P>=0.5, tie P>=0.25), which demotes small-field races whose
extreme point v_j is noise (winner's curse). min@3:00 =
180*(exp(v_j)-1) = minutes vs the average race for a 3:00:00 runner.

slice: WA_M_14-25_mrc2; min_n = 0; n_top = 20; model = full_nu8p00

### WA_M_14-25_mrc2 -- slowest 20 (headline P>=0.5, tie P>=0.25)

| rank | race | country | n_j | v_j | v 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|
| 1 | nyc_marathon 2022 | USA | 102 | +0.0586 | [+0.0476, +0.0681] | +10.9 | 1.00 | [1, 3] | headline |
| 2 | stockholm_marathon 2024 | SWE | 32 | +0.0536 | [+0.0459, +0.0619] | +9.9 | 1.00 | [1, 3] | headline |
| 3 | boston_marathon 2018 | USA | 60 | +0.0481 | [+0.0383, +0.0568] | +8.9 | 1.00 | [1, 5] | headline |
| 4 | boston_marathon 2016 | USA | 54 | +0.0371 | [+0.0308, +0.0466] | +6.8 | 1.00 | [3, 7] | headline |
| 5 | london_marathon 2018 | GBR | 225 | +0.0360 | [+0.0301, +0.0418] | +6.6 | 1.00 | [4, 9] | headline |
| 6 | sydney_marathon 2023 | AUS | 20 | +0.0310 | [+0.0117, +0.0461] | +5.7 | 0.94 | [3, 37] | headline |
| 7 | boston_marathon 2017 | USA | 95 | +0.0292 | [+0.0214, +0.0370] | +5.3 | 1.00 | [5, 14] | headline |
| 8 | nyc_marathon 2017 | USA | 95 | +0.0258 | [+0.0188, +0.0316] | +4.7 | 0.99 | [6, 18] | headline |
| 9 | nyc_marathon 2023 | USA | 99 | +0.0252 | [+0.0185, +0.0318] | +4.6 | 0.97 | [7, 21] | headline |
| 10 | nyc_marathon 2016 | USA | 81 | +0.0244 | [+0.0183, +0.0327] | +4.4 | 0.97 | [6, 21] | headline |
| 11 | stockholm_marathon 2023 | SWE | 38 | +0.0229 | [+0.0160, +0.0327] | +4.2 | 0.90 | [6, 25] | headline |
| 12 | nyc_marathon 2015 | USA | 65 | +0.0223 | [+0.0161, +0.0302] | +4.1 | 0.79 | [8, 29] | headline |
| 13 | nyc_marathon 2025 | USA | 118 | +0.0223 | [+0.0155, +0.0291] | +4.1 | 0.89 | [7, 30] | headline |
| 14 | nyc_marathon 2014 | USA | 66 | +0.0212 | [+0.0122, +0.0302] | +3.9 | 0.68 | [7, 37] | headline |
| 15 | stockholm_marathon 2022 | SWE | 31 | +0.0207 | [+0.0128, +0.0306] | +3.8 | 0.67 | [7, 35] | headline |
| 16 | prague_marathon 2024 | CZE | 25 | +0.0202 | [+0.0101, +0.0299] | +3.7 | 0.59 | [8, 43] | headline |
| 17 | nyc_marathon 2024 | USA | 134 | +0.0202 | [+0.0152, +0.0263] | +3.7 | 0.72 | [11, 28] | headline |
| 18 | dublin_marathon 2017 | IRL | 29 | +0.0192 | [+0.0121, +0.0299] | +3.5 | 0.49 | [7, 36] | tie |
| 19 | copenhagen_marathon 2018 | DEN | 20 | +0.0183 | [-0.0014, +0.0384] | +3.3 | 0.54 | [4, 95] | headline |
| 20 | edinburgh_marathon 2022 | GBR | 27 | +0.0169 | [-0.0040, +0.0452] | +3.1 | 0.54 | [4, 111] | headline |
| 21 | stockholm_marathon 2025 | SWE | 43 | +0.0169 | [+0.0121, +0.0255] | +3.1 | 0.39 | [11, 35] | tie |
| 22 | milan_marathon 2023 | ITA | 29 | +0.0168 | [+0.0026, +0.0299] | +3.1 | 0.31 | [7, 71] | tie |
| 24 | chicago_marathon 2021 | USA | 66 | +0.0166 | [+0.0097, +0.0260] | +3.0 | 0.28 | [9, 43] | tie |
| 26 | boston_marathon 2021 | USA | 109 | +0.0165 | [+0.0094, +0.0226] | +3.0 | 0.29 | [14, 44] | tie |

### WA_M_14-25_mrc2 -- fastest 20 (headline P>=0.5, tie P>=0.25)

| rank | race | country | n_j | v_j | v 95% CI | min@3:00 | P(topN) | rank 95% | tier |
|---|---|---|---|---|---|---|---|---|---|
| 1 | tokyo_marathon 2015 | JPN | 69 | -0.0280 | [-0.0350, -0.0197] | -5.0 | 0.99 | [1, 9] | headline |
| 2 | zurich_marathon 2025 | SUI | 20 | -0.0227 | [-0.0312, -0.0077] | -4.0 | 0.82 | [1, 65] | headline |
| 3 | valencia_marathon 2023 | ESP | 659 | -0.0212 | [-0.0238, -0.0190] | -3.8 | 1.00 | [2, 12] | headline |
| 4 | seville_marathon 2020 | ESP | 99 | -0.0209 | [-0.0253, -0.0145] | -3.7 | 0.95 | [2, 30] | headline |
| 5 | valencia_marathon 2020 | ESP | 49 | -0.0205 | [-0.0271, -0.0161] | -3.7 | 0.96 | [2, 22] | headline |
| 6 | berlin_marathon 2014 | GER | 73 | -0.0199 | [-0.0264, -0.0135] | -3.5 | 0.78 | [2, 36] | headline |
| 7 | tokyo_marathon 2018 | JPN | 114 | -0.0198 | [-0.0253, -0.0154] | -3.5 | 0.94 | [2, 26] | headline |
| 8 | hannover_marathon 2025 | GER | 40 | -0.0198 | [-0.0287, -0.0135] | -3.5 | 0.86 | [1, 33] | headline |
| 9 | seville_marathon 2022 | ESP | 130 | -0.0194 | [-0.0236, -0.0147] | -3.5 | 0.90 | [2, 28] | headline |
| 10 | hannover_marathon 2022 | GER | 36 | -0.0191 | [-0.0276, -0.0124] | -3.4 | 0.77 | [2, 43] | headline |
| 11 | chicago_marathon 2014 | USA | 78 | -0.0182 | [-0.0238, -0.0117] | -3.3 | 0.73 | [2, 43] | headline |
| 12 | dublin_marathon 2025 | IRL | 55 | -0.0181 | [-0.0266, -0.0098] | -3.2 | 0.65 | [2, 54] | headline |
| 13 | seville_marathon 2024 | ESP | 297 | -0.0180 | [-0.0223, -0.0147] | -3.2 | 0.86 | [4, 27] | headline |
| 14 | amsterdam_marathon 2025 | NED | 114 | -0.0175 | [-0.0239, -0.0132] | -3.1 | 0.61 | [3, 43] | headline |
| 15 | berlin_marathon 2015 | GER | 164 | -0.0174 | [-0.0214, -0.0140] | -3.1 | 0.72 | [6, 33] | headline |
| 16 | valencia_marathon 2019 | ESP | 97 | -0.0166 | [-0.0214, -0.0106] | -3.0 | 0.48 | [5, 53] | tie |
| 17 | frankfurt_marathon 2015 | GER | 56 | -0.0164 | [-0.0225, -0.0092] | -2.9 | 0.39 | [5, 62] | tie |
| 18 | eindhoven_marathon 2022 | NED | 49 | -0.0161 | [-0.0235, -0.0092] | -2.9 | 0.61 | [2, 57] | headline |
| 19 | valencia_marathon 2025 | ESP | 443 | -0.0161 | [-0.0188, -0.0122] | -2.9 | 0.31 | [12, 42] | tie |
| 20 | tokyo_marathon 2022 | JPN | 127 | -0.0159 | [-0.0215, -0.0103] | -2.8 | 0.42 | [4, 50] | tie |
| 21 | amsterdam_marathon 2021 | NED | 79 | -0.0159 | [-0.0231, -0.0110] | -2.8 | 0.38 | [5, 49] | tie |
| 22 | amsterdam_marathon 2023 | NED | 113 | -0.0155 | [-0.0195, -0.0112] | -2.8 | 0.37 | [10, 46] | tie |
| 23 | berlin_marathon 2019 | GER | 249 | -0.0152 | [-0.0185, -0.0124] | -2.7 | 0.25 | [13, 41] | tie |
| 24 | chicago_marathon 2019 | USA | 147 | -0.0152 | [-0.0203, -0.0104] | -2.7 | 0.39 | [7, 55] | tie |
| 27 | boston_marathon 2014 | USA | 87 | -0.0141 | [-0.0201, -0.0062] | -2.5 | 0.33 | [10, 74] | tie |
| 33 | yorkshire_marathon 2019 | GBR | 26 | -0.0131 | [-0.0281, -0.0015] | -2.3 | 0.35 | [1, 102] | tie |
| 46 | paris_marathon 2015 | FRA | 20 | -0.0109 | [-0.0305, +0.0141] | -2.0 | 0.34 | [1, 166] | tie |
