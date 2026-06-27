# Cross-slice big movers (gauge-free, beyond-noise series rank shifts)

r_j = (v_A - v_B) with {1, t, t^2} removed over the shared races =
the part of the difference no gauge convention can produce; r > 0 =
race relatively SLOWER in A than B. Series contrast = median of edition
r_j. Flagged: |180*r| >= 0.5 min at 3:00 AND z = |r|/sd_boot >= 2.0.
P(top-N) is deliberately NOT used (N-dependent, jittery).

| comparison | series | ctry | k | contrast r | min@3:00 | z | rank A->B (d) |
|---|---|---|---|---|---|---|---|
| MvsW_ALL_14-25_mrc2 | berlin_marathon | GER | 11 | -0.0051 | -0.91 | 6.3 | 14->21 (+7) |
| MvsW_ALL_14-25_mrc2 | boston_marathon | USA | 11 | +0.0073 | +1.31 | 5.8 | 36->34 (-2) |
| MvsW_ALL_14-25_mrc2 | chicago_marathon | USA | 11 | -0.0041 | -0.73 | 4.5 | 21->24 (+3) |
| MvsW_ALL_14-25_mrc2 | valencia_marathon | ESP | 9 | -0.0043 | -0.78 | 3.3 | 8->15 (+7) |
| ALLvsPo10_W_14-25_mrc2 | boston_marathon | USA | 10 | -0.0071 | -1.27 | 3.1 | 22->25 (+3) |
| MvsW_ALL_14-25_mrc2 | hannover_marathon | GER | 6 | -0.0076 | -1.37 | 2.9 | 4->19 (+15) |
| ALLvsPo10_W_14-25_mrc2 | paris_marathon | FRA | 10 | +0.0115 | +2.08 | 2.8 | 14->10 (-4) |
| MvsW_ALL_14-25_mrc2 | milan_marathon | ITA | 4 | -0.0055 | -0.99 | 2.6 | 33->36 (+3) |
| ALLvsPo10_M_14-25_mrc2 | stockholm_marathon | SWE | 8 | -0.0122 | -2.19 | 2.6 | 19->28 (+9) |
| MvsW_Po10_14-25_mrc2 | chicago_marathon | USA | 10 | -0.0060 | -1.08 | 2.6 | 13->16 (+3) |
| MvsW_ALL_14-25_mrc2 | prague_marathon | CZE | 10 | -0.0107 | -1.92 | 2.6 | 23->35 (+12) |
| MvsW_ALL_14-25_mrc2 | seville_marathon | ESP | 7 | -0.0046 | -0.82 | 2.5 | 3->4 (+1) |
| MvsW_Po10_14-25_mrc2 | berlin_marathon | GER | 11 | -0.0043 | -0.77 | 2.5 | 11->15 (+4) |
| MvsW_ALL_14-25_mrc2 | malaga_marathon | ESP | 7 | +0.0061 | +1.11 | 2.4 | 20->3 (-17) |
| MvsW_ALL_14-25_mrc2 | rome_marathon | ITA | 4 | -0.0052 | -0.93 | 2.3 | 31->37 (+6) |
| MvsW_ALL_14-25_mrc2 | cape_town_marathon | RSA | 3 | +0.0217 | +3.90 | 2.2 | 38->22 (-16) |
| ALLvsPo10_W_14-25_mrc2 | nyc_marathon | USA | 11 | -0.0045 | -0.81 | 2.2 | 25->24 (-1) |
| MvsW_ALL_14-25_mrc2 | cologne_marathon | GER | 6 | +0.0047 | +0.85 | 2.1 | 25->8 (-17) |
| MvsW_ALL_14-25_mrc2 | melbourne_marathon | AUS | 11 | +0.0090 | +1.61 | 2.1 | 40->28 (-12) |
| ALLvsPo10_M_14-25_mrc2 | paris_marathon | FRA | 11 | +0.0051 | +0.91 | 2.1 | 10->11 (+1) |
| MvsW_ALL_14-25_mrc2 | oslo_marathon | NOR | 8 | +0.0070 | +1.25 | 2.1 | 34->27 (-7) |
| ALLvsPo10_M_14-25_mrc2 | munich_marathon | GER | 3 | -0.0204 | -3.67 | 2.1 | 16->30 (+14) |
