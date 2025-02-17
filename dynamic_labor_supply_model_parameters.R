## 1. time and periods ----
n_periods_year <- 4
t1_age <- 29 # at least 1 year younger than treatment age
bi_age <- 30
retirement_age <- 65
life_expectancy <- 80
n_bi_years <- 3

terminal_period <- (life_expectancy - t1_age) * n_periods_year
retirement_period <- (retirement_age - t1_age - 1) * n_periods_year + 1
bi_start_period <- (bi_age - t1_age) * n_periods_year + 1
bi_end_period <- bi_start_period + 3*n_periods_year -1
bl_period <- round(bi_start_period - 0.51*n_periods_year)
ml_period <- round(bi_start_period + 1.51*n_periods_year)
el_period <- round(bi_start_period + 2.51*n_periods_year)

## 2. asset space and spline function ----
A_min <- -1e5
A_max <- 1e6
A_nodes_curv <- 0.3
n_nodes <- 50

A_nodes = seq(-(abs(A_min)^A_nodes_curv), A_max^A_nodes_curv, length.out = n_nodes)
A_nodes[A_nodes>0] = A_nodes[A_nodes>0] ^ (1/A_nodes_curv)
A_nodes[A_nodes<0] = -(abs(A_nodes[A_nodes<0]) ^ (1/A_nodes_curv))
A_nodes[1] = A_min # this is a technical necessity because of rounding
A_nodes[n_nodes] = A_max # this is a technical necessity because of rounding


## 3. rolling windows ----
W_year_min = 0
W_year_max = retirement_age - t1_age -1
W_curve = 0.5
W_length = 30
Ws_yearly = (seq(W_year_min^W_curve, W_year_max^W_curve, length.out = W_length)^(1/W_curve))

Ws = unique(round(Ws_yearly*n_periods_year))
Ws <- subset(Ws, Ws != 0)


## 4. employment model ----
π0 = matrix(c(0.6, 0.4))
unemp_to_emp_half_yearly = 1/3
unemp_to_emp_optimistic_half_yearly = 0.7

Π = solve_transition_matrix(π0, unemp_to_emp_half_yearly, 4)
Π_optimistic_whenunemp = solve_transition_matrix_optimistic(
  unemp_to_emp_optimistic_half_yearly,
  unemp_to_emp_half_yearly,
  n_periods_year
)



