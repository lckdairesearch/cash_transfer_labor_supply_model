# Treatment Effects and SE ----
c_effect = tribble(
  ~"source", ~"month", ~"beta", ~"std.error", ~"controlmean",
  "src", 19, 260.28775, 58.415652, 4222.65,
  "src", 30, 217.41985, 62.890364, 4349.464,
  "qualtrics", 4, 184.01618, 54.616082, 3632.286,
  "qualtrics", 8, 269.86747, 55.13548, 3686.903,
  "qualtrics", 12, 291.44042, 56.545459, 3661.506,
  "qualtrics", 15, 325.53306, 64.546238, 3887.569,
  "qualtrics", 28, 417.31903, 70.385203, 3959.216,
  "qualtrics", 31, 332.27861, 60.28486, 3865.902,
  "qualtrics", 36, 387.24385, 68.311982, 3848.909,
  "qualtrics", 38, 237.49192, 73.639675, 4116.65,
  "qualtrics", 43, 295.05154, 80.906357, 4157.842,
)

c_effect_paper = tribble(
  ~"src", ~"month", ~"beta", ~"std.error", ~"controlmean",
  "qualtrics", 6, 279.6732, 64.51266, 3660.204,
  "qualtrics", 18, 323.5165, 71.15285, 4068.961,
  "qualtrics", 30, 356.1528, 66.45753, 4050.373,
)


l_effect = tribble(
  ~"source", ~"month", ~"beta", ~"std.error", ~"controlmean",
  "qualtrics", 4, -0.07619937, 0.6729077, 24.66015,
  "qualtrics", 7, -0.50849787, 0.6860296, 25.1977,
  "qualtrics", 10, -0.52675998, 0.6977845, 25.56048,
  "qualtrics", 13, -1.02299714, 0.6949056, 27.12273,
  "qualtrics", 16, -0.82060034, 0.711247, 27.84197,
  "qualtrics", 19, -0.62256158, 0.6597204, 28.22114,
  "qualtrics", 22, -2.12156512, 0.7536094, 28.47107,
  "qualtrics", 25, -1.90201086, 0.7667052, 29.63133,
  "qualtrics", 28, -1.71868631, 0.7735355, 30.18108,
  "qualtrics", 31, -1.59477216, 0.7758373, 30.129,
  "qualtrics", 34, -0.80062791, 0.8001348, 27.24345,
  "qualtrics", 37, -0.97429855, 0.7880694, 26.65817,
  "qualtrics", 40, -0.41032125, 0.8344483, 26.51707,
)

l_effect_paper = tribble(
  ~"source", ~"month", ~"beta", ~"std.error", ~"controlmean",
  "src", 19, -0.786582, 0.710166, 28.71634,
  "src", 30, -1.554392, 0.771105, 31.13073,
)


A_effect = tribble(
  ~"type", ~"month", ~"value", ~"std.error", ~"controlmean",
  "control_raw", 18, -5351.629, 67639.52, NA,
  "control_raw", 30, -2241.443, 75045.11, NA,
  "treat_effect", 18, -651, 2740, NA,
  "treat_effect", 30, -2680, 2790, NA 
)


income_effect_admin <- tribble( 
  ~"source", ~"month", ~"beta", ~"std.error", ~"controlmean", 
  "admin", 3, -313.092, 219.3046, 3507.515, 
  "admin", 6, -383.6166, 236.1418, 3850.149, 
  "admin", 9, -118.3641, 251.7471, 4008.626, 
  "admin", 12, -187.8067, 284.2417, 4684.816, 
  "admin", 15, -253.1292, 274.5954, 4552.06, 
  "admin", 18, -370.8209, 294.3916, 4986.341, 
  "admin", 21, -429.2238, 315.7866, 5442.901, 
  "admin", 24, -631.0487, 325.4435, 5668.161, 
  "admin", 27, -662.7, 326.8954, 5625.651, 
  "admin", 30, -577.2886, 329.4384, 5719.71, 
  "admin", 33, -511.5341, 330.5609, 5762.391, 
  "admin", 36, -503.8665, 329.6042, 5888.784, 
  "admin", 39, -422.227, 332.5881, 5855.033 
)



model_hh_size = 1.7
c_effect_underreport_factor =  (0.33/0.74)
c_effect[c(3,4)] <- c_effect[c(3,4)]/c_effect_underreport_factor
c_effect[-c(1, 2)] <- c_effect[-c(1, 2)]*12/n_periods_year
c_effect_paper[c(3,4)] <- c_effect_paper[c(3,4)]/c_effect_underreport_factor
c_effect_paper[-c(1, 2)] <- c_effect_paper[-c(1, 2)]*12/n_periods_year
l_effect[-c(1, 2)] <- l_effect[-c(1, 2)]*52/n_periods_year*model_hh_size
l_effect_paper[-c(1, 2)] <- l_effect_paper[-c(1, 2)]*52/n_periods_year*model_hh_size
income_effect_admin[-c(1, 2)] <- income_effect_admin[-c(1, 2)]*12/n_periods_year*model_hh_size
