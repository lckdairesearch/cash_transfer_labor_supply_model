# 🌝 PARAMETERS ----
compute_B <- function(
    σ, η, c, l, w, γc
){
  B <- (c-γc)^(-σ) * l^(-(1/η)) * w
  return(B)
}

compute_period_B_from_yearly <- function(
    B_yearly, σ, η,
    .n_periods_year = n_periods_year
){
  B <- ((1/.n_periods_year)^(-(σ+(1/η))))*B_yearly
  return(B)
}

get_B_grid = function(
    σ, η, n_Bs_out,
    c_bl, l_bl
){
  c_bl_grid <- seq(c_bl[1] - 0.5*c_bl[2], c_bl[1] + 0.5*c_bl[2], length.out = 5)
  l_bl_grid <- seq(l_bl[1] - 0.5*l_bl[2], l_bl[1] + 0.5*l_bl[2], length.out = 5)
  c_l_bl_grid <- expand.grid(c_bl = c_bl_grid, l_bl = l_bl_grid)
  B_grid <- winsorize(unlist(map2(c_l_bl_grid$c_bl, c_l_bl_grid$l_bl, \(c_bl, l_bl) compute_B(σ, η, c_bl, l_bl, w_bl, γc))))
  B_grid <- seq(min(B_grid), max(B_grid), length.out = (n_Bs_out - 1))
  B <- compute_B(σ, η, c_bl[1], l_bl[1], w_bl, γc)
  B_grid <- sort(unique(append(B_grid,B)))
  
  return(B_grid)
}

solve_transition_matrix = function(π0, unemp_to_emp_half_yearly, n_periods_year) {
  
  emp_to_emp <- function(unemp_to_emp) {
    return((π0[1] - π0[2]*unemp_to_emp)/π0[1])
  }
  
  Π = function(unemp_to_emp) {
    return(matrix(c(
      emp_to_emp(unemp_to_emp), (1-emp_to_emp(unemp_to_emp)),
      unemp_to_emp, (1-unemp_to_emp)
    ), nrow = 2, byrow =T))
  }
  
  Π_half_yearly = Π(unemp_to_emp_half_yearly)
  
  Π = expm(logm(Π_half_yearly) / (n_periods_year/2))
  
  return(Π)
  
}

solve_transition_matrix_optimistic = function(unemp_to_emp_optimistic_half_yearly, unemp_to_emp_half_yearly, n_periods_year){
  
  # here we assume that they have correct belief about their employed transition probability at the half year mark and nothing more
  
  emp_to_emp <- function(unemp_to_emp) {
    return((π0[1] - π0[2]*unemp_to_emp)/π0[1])
  }
  
  Π_optimistic_half_yearly = matrix(c(
    emp_to_emp(unemp_to_emp_half_yearly), (1-emp_to_emp(unemp_to_emp_half_yearly)),
    unemp_to_emp_optimistic_half_yearly, (1-unemp_to_emp_optimistic_half_yearly)
  ), nrow = 2, byrow =T)
  
  Π_optimistic = expm(logm(Π_optimistic_half_yearly) / (n_periods_year/2))
  
  
  return(Π_optimistic)
  
}

define_TRt <- function(
    t, treat_status, 
    treat_transfer_monthly, control_transfer_monthly,
    .bi_start_period = bi_start_period, .bi_end_period = bi_end_period, 
    .n_periods_year = n_periods_year
) {
  if (t %in% seq(.bi_start_period, .bi_end_period) & treat_status == 'treat') {
    TRt = treat_transfer_monthly * 12/.n_periods_year
  } else if (t %in% seq(.bi_start_period, .bi_end_period) & treat_status == 'control') {
    TRt = control_transfer_monthly * 12 /.n_periods_year
  } else {
    TRt = 0
  }
  return(TRt)
}




# 🌊 UTILITY FUNCTIONS ----
compute_u <- function(
    c, σ,
    .γc = γc
){
  c_net = max(c-.γc, 1e-5) # c_net must be non-negative
  utility <- ifelse(σ == 1,
    {log(c_net)},
    {(1-σ)^(-1)*(c_net)^(1-σ)}
  )
  return(utility)
}

compute_v <- function(
    l, η
){
  disutility <- (η/(1+η))*l^((1+η)/η)
  return(disutility)
}

compute_Vt <- function(
    ct, lt, Vnextt,
    σ, η, B, δ,
    qhyp = 0, β = 1,
    .compute_u = compute_u, .compute_v = compute_v
){
  Vt <- ifelse(qhyp == 0,
    {.compute_u(ct, σ) - B*.compute_v(lt, η) + δ*Vnextt},
    {.compute_u(ct, σ) - B*.compute_v(lt, η) + β*δ*Vnextt}
  )
  return(Vt)
}

# 🧠 CLOSED FORM SOLUTION FOR ct lt ----
compute_cstart <- function(
    lstart, wt,
    σ, η, B,
    .γc = γc
){
  cstart <- (B*lstart^(1/η)/(wt))^(-1/σ) + .γc
  return(cstart)
}



solve_cstart_lstart <-function(
    t, At, Anextt, It, TRt, wt, lmax,
    σ, η, B, 
    .retirement_period = retirement_period, .r = r, .r_borrowing = r_borrowing,
    .compute_cstart = compute_cstart
){
  if (t < .retirement_period) {
    # if Working
    budget_constraint <- function(lstart) {
      (At - .compute_cstart(lstart, wt, σ, η, B) + wt*lstart + It + TRt) - Anextt/(1+.r)*(Anextt >= 0) - Anextt*(1+.r_borrowing)*(Anextt < 0)
    }
    c(lstart, cstart) %<-% tryCatch({
      # working case 1. interior solution
      lstart <- (uniroot(budget_constraint, interval = c(0, lmax)))$root
      cstart <- .compute_cstart(lstart, wt, σ, η, B)
      c(lstart, cstart) %<-% list(lstart, cstart)
    }, 
    error = function(e) {
      # working case 2. edge cases
      if (budget_constraint(lmax) < 0) {
        # working case 2a. consumption from FOC conditions at lmax is not affordable
        lstart <- lmax
        cstart <- At + It + TRt + wt*lstart - Anextt/(1+.r)*(Anextt >= 0) - Anextt*(1+.r_borrowing)*(Anextt < 0)
      } else if (budget_constraint(lmax) > 0) {
        # working case 2b. consumption from FOC conditions at lt = 0 still does not exhaust budget
        lstart <- 0
        cstart <- At + It + TRt + wt*lstart - Anextt/(1+.r)*(Anextt >= 0) - Anextt*(1+.r_borrowing)*(Anextt < 0)
      }
      return(list(lstart, cstart))
    })
  } else{
    
    # else Retired
    lstart <- 0
    cstart <- At + It + TRt - Anextt/(1+.r)*(Anextt >= 0) - Anextt*(1+.r_borrowing)*(Anextt < 0)
  }
  
  return(c(cstart = cstart, lstart = lstart))
}

# 🧠 Vtmax at At ----
Vt_atAt_givenAnextt <- function(
    t, At, Anextt, models,
    It, TRt, wt, lmax,
    σ, η, B, δ,
    qhyp = 0, β = 1,
    .A_max = A_max, .A_min = A_min, .n_nodes = n_nodes,
    .solve_cstart_lstart = solve_cstart_lstart, .compute_Vt = compute_Vt
){
  c(ct, lt) %<-% .solve_cstart_lstart(t, At, Anextt, It, TRt, wt, lmax, σ, η, B)
  Vnextt <- predict(models[[t+1]], Anextt)$y
  Vt <- .compute_Vt(ct, lt, Vnextt, σ, η, B, δ, qhyp, β)
  return(Vt)
}

Vtmax_atAt <- function(
    t, At, models,
    It, TRt, wt, lmax,
    σ, η, B, δ,
    qhyp = 0, β = 1,
    .r = r, .r_borrowing = r_borrowing, .retirement_period = retirement_period, 
    .A_max = A_max, .A_min = A_min, .n_nodes = n_nodes,
    .Vt_atAt_givenAnextt = Vt_atAt_givenAnextt, .compute_cstart = compute_cstart
){
  Vt_atAt_givenAnextt_partial <- partial(
    .Vt_atAt_givenAnextt,
    t = t, At = At, models = models,
    It = It, TRt = TRt, wt = wt, lmax = lmax,
    σ = σ, η = η, B = B, δ = δ,
    qhyp = qhyp, β = β,
  )
  lb <- models$fit$min
  lb <- max(lb, .A_min)
  ub <- At + It + TRt + wt*lmax*(t < .retirement_period)
  ub <- ub * ((1+.r)*(ub>=0) + (1+.r_borrowing)*(ub<0))
  ub <- min(ub, .A_max)
  if (ub > A_min & lb < A_max) {
    # normal case
    optimized =  optimize(Vt_atAt_givenAnextt_partial, interval = c(.A_min, ub), maximum = TRUE)
    Vtmax = optimized$objective
    Anextt = optimized$maximum
    
  } else if (ub <= A_min) {
    # edge case 
    Vtmax = Vt_atAt_givenAnextt_partial(A_min)
    Anextt = .A_min
  } else if (lb >= A_max) {
    # edge case
    Vtmax = Vt_atAt_givenAnextt_partial(A_max)
    Anextt = .A_max
  }
  return(c(Vtmax = Vtmax, Anextt = Anextt))
}



# 🌵 STANDARD MODEL ----
## Solve ----
solve_dynamic_labor_supply_model <- function(
    last_period,
    It, w_bl, lmax,
    σ, η, B, δ,
    treat_transfer_monthly, control_transfer_monthly,
    .ω = ω,.bl_period = bl_period, 
    .γc = γc, .r = r, .r_borrowing = r_borrowing,
    .bi_start_period = bi_start_period, .bi_end_period = bi_end_period, .retirement_period = retirement_period,
    .n_nodes = n_nodes, .A_nodes = A_nodes, .A_min = A_min,
    .define_TRt = define_TRt,.compute_Vt = compute_Vt, .Vtmax_atAt = Vtmax_atAt,
    .solve_cstart_lstart = solve_cstart_lstart
){
  
  
  # 1. Not in Experiment
  models_normal <- list()
  # 1.1 last period
  t <- last_period
  wt <- w_bl*(1+.ω)^(t - .bl_period)
  TRt <- .define_TRt(t, 'normal', treat_transfer_monthly, control_transfer_monthly)
  bundles_nodes <- map(.A_nodes, 
                       \(At) 
                       .solve_cstart_lstart(
                         t, At, Anextt = 0, It, TRt, wt, lmax,
                         σ, η, B))
  Vt_nodes <- map2(
    map_dbl(bundles_nodes, 'cstart'),
    map_dbl(bundles_nodes, 'lstart'),
    \(ct, lt)
    .compute_Vt(ct, lt, Vnextt = 0, σ, η, B, δ))
  At_min <- -It + .γc
  models_normal[[t]] <- ss(.A_nodes[.A_nodes > At_min], Vt_nodes[.A_nodes > At_min], method="OCV", all.knots = T, xmin = At_min)
  models_normal[[t+1]] <- ss(.A_nodes, rep(0, length(.A_nodes)), method="OCV", all.knots = T)
  
  # 1.2 Period from last_period-1 to 1
  for (t in seq(last_period-1, 1)) {
    wt <- w_bl*(1+.ω)^(t - .bl_period)
    TRt <- .define_TRt(t, 'normal', treat_transfer_monthly, control_transfer_monthly)
    Vtmax_atAt_partial <- partial(
      .Vtmax_atAt,
      t = t, models = models_normal,
      It = It, TRt = TRt, wt = wt, lmax = lmax,
      σ=σ, η =η, B = B, δ = δ)
    Vt_nodes <- map_dbl(map(.A_nodes, Vtmax_atAt_partial), 'Vtmax')
    Anextt_min <- models_normal[[t+1]]$fit$min
    At_min <- Anextt_min/((1+.r)*(Anextt_min >=0) + (1+.r_borrowing)*(Anextt_min <0)) - It - TRt - wt*lmax*(t < .retirement_period)  + .γc
    At_min <- max(At_min, .A_min)
    models_normal[[t]] <- ss(.A_nodes[.A_nodes > At_min], Vt_nodes[.A_nodes > At_min], method="OCV", all.knots = T, xmin = At_min)
  }
  
  # 2. Control Group 
  models_control <- models_normal
  # From bi_end_period to bi_start_period
  for (t in seq(min(last_period, .bi_end_period), min(last_period, .bi_start_period))) {
    wt <- w_bl*(1+.ω)^(t - .bl_period)
    TRt <- .define_TRt(t, 'control', treat_transfer_monthly, control_transfer_monthly)
    Vtmax_atAt_partial <- partial(
      .Vtmax_atAt,
      t = t, models = models_control,
      It = It, TRt = TRt, wt = wt, lmax = lmax,
      σ=σ, η =η, B = B, δ = δ)
    Vt_nodes <- map_dbl(map(.A_nodes, Vtmax_atAt_partial), 'Vtmax')
    Anextt_min <- models_control[[t+1]]$fit$min
    At_min <- Anextt_min/((1+.r)*(Anextt_min >=0) + (1+.r_borrowing)*(Anextt_min <0)) - It - TRt - wt*lmax*(t < .retirement_period)  + .γc
    At_min <- max(At_min, .A_min)
    models_control[[t]] <- ss(.A_nodes[.A_nodes > At_min], Vt_nodes[.A_nodes > At_min], method="OCV", all.knots = T, xmin = At_min)
  }
  
  # 3. Treatment Group 
  models_treat <- models_normal
  # From bi_end_period to bi_start_period
  for (t in seq(min(last_period, .bi_end_period), min(last_period, .bi_start_period))) {
    wt <- w_bl*(1+.ω)^(t - .bl_period)
    TRt = .define_TRt(t, 'treat', treat_transfer_monthly, control_transfer_monthly)
    Vtmax_atAt_partial <- partial(
      .Vtmax_atAt,
      t = t, models = models_normal,
      It = It, TRt = TRt, wt = wt, lmax = lmax,
      σ=σ, η =η, B = B, δ = δ)
    Vt_nodes <- map_dbl(map(.A_nodes, Vtmax_atAt_partial), 'Vtmax')
    Anextt_min <- models_treat[[t+1]]$fit$min
    At_min <- Anextt_min/((1+.r)*(Anextt_min >=0) + (1+.r_borrowing)*(Anextt_min <0)) - It - TRt - wt*lmax*(t < .retirement_period)  + .γc
    At_min <- max(At_min, .A_min)
    models_treat[[t]] <- ss(.A_nodes[.A_nodes > At_min], Vt_nodes[.A_nodes > At_min], method="OCV", all.knots = T, xmin = At_min)
  }
  
  return(list(normal = models_normal, control = models_control, treat = models_treat))
}


## Paths ----
get_paths <- function(
    A1, n_periods, models_bystatus,
    w_bl, It, lmax,
    σ, η, B, δ, 
    treat_transfer_monthly, control_transfer_monthly,
    qhyp = 0, β = 1,
    .ω = ω,
    .bl_period = bl_period,
    .bi_start_period = bi_start_period, .bi_end_period = bi_end_period,
    .Vtmax_atAt = Vtmax_atAt, .solve_cstart_lstart = solve_cstart_lstart,
    .define_TRt = define_TRt
){
  # Variables
  c_path_normal = list()
  c_path_control = list()
  c_path_treat = list()
  l_path_normal = list()
  l_path_control = list()
  l_path_treat = list()
  A_path_normal = list()
  A_path_control = list()
  A_path_treat = list()
  
  
  find_controlvar_att = function (t, At, models, It, TRt, wt, lmax, σ, η, B, δ, qhyp, β) {
    ## find Anextt
    c(Vt, Anextt) %<-% .Vtmax_atAt(t, At, models, It, TRt, wt, lmax, σ, η, B, δ, qhyp, β)
    ## solve optimal bundle
    bundle = .solve_cstart_lstart(t, At, Anextt, It, TRt, wt, lmax, σ, η, B)
    cstart = bundle['cstart']
    lstart = bundle['lstart']
    return(c(lstart, cstart, Anextt))
  }
  
  # Initial asset   
  At_normal = A1
  At_control = A1
  At_treat = A1
  
  # Paths
  ## Normal
  for (t in seq(1, n_periods)) {
    wt <- w_bl*(1+.ω)^(t - .bl_period)
    TRt <- .define_TRt(t, 'normal', treat_transfer_monthly, control_transfer_monthly)
    c(lt_normal, ct_normal, Anextt_normal) %<-%  find_controlvar_att(t, At_normal, models_bystatus$normal, It, TRt, wt, lmax, σ, η, B, δ, qhyp, β)
    l_path_normal[t] <- lt_normal
    c_path_normal[t] <- ct_normal
    A_path_normal[t] <- At_normal
    At_normal <- Anextt_normal
  }
  
  
  ## Control
  for (t in seq(1, n_periods)) {
    wt <- w_bl*(1+.ω)^(t - .bl_period)
    TRt <- .define_TRt(t, 'control', treat_transfer_monthly, control_transfer_monthly)
    ### Should expect no change in TR before start of benefit
    if(t < .bi_start_period) models_control_att <- models_bystatus$normal
    else models_control_att <- models_bystatus$control
    c(lt_control, ct_control, Anextt_control) %<-%  find_controlvar_att(t, At_control, models_control_att, It, TRt, wt, lmax, σ, η, B, δ, qhyp, β)
    l_path_control[t] <- lt_control
    c_path_control[t] <- ct_control
    A_path_control[t] <- At_control
    At_control <- Anextt_control
  }
  
  ## Treat
  for (t in seq(1, n_periods)) {
    wt <- w_bl*(1+.ω)^(t - .bl_period)
    TRt <- .define_TRt(t, 'treat', treat_transfer_monthly, control_transfer_monthly)
    ### Should expect no change in TR before start of benefit
    if(t < .bi_start_period) models_treat_att <- models_bystatus$normal
    else models_treat_att <- models_bystatus$treat
    c(lt_treat, ct_treat, Anextt_treat) %<-%  find_controlvar_att(t, At_treat, models_treat_att, It, TRt, wt, lmax, σ, η, B, δ, qhyp, β)
    l_path_treat[t] <- lt_treat
    c_path_treat[t] <- ct_treat
    A_path_treat[t] <- At_treat
    At_treat <- Anextt_treat
  }
  
  return(list(c_path_normal, c_path_control, c_path_treat, l_path_normal, l_path_control, l_path_treat, A_path_normal, A_path_control, A_path_treat))
}

# 🪟 ROLLING WINDOWS MODEL ----
## Solve ----
solve_dynamic_labor_supply_rolling_windows_model <- function(
    last_period,
    It, w_bl, lmax,
    σ, η, B, δ,
    treat_transfer_monthly, control_transfer_monthly,
    .solve_dynamic_labor_supply_model = solve_dynamic_labor_supply_model
){
  periods = seq(1, last_period , 1)
  coeffs_windows = map(periods, \(period) .solve_dynamic_labor_supply_model(period, It, w_bl, lmax, σ, η, B, δ, treat_transfer_monthly, control_transfer_monthly))
  return(coeffs_windows)
}

get_paths_rolling_windows <- function(
    A1, n_periods, coeffs_windows,
    w_bl, It, lmax,
    σ, η, B, δ, W,
    treat_transfer_monthly, control_transfer_monthly,
    qhyp = 0, β = 1,
    .ω = ω,
    .bl_period = bl_period,
    .bi_start_period = bi_start_period, .bi_end_period = bi_end_period,
    .Vtmax_atAt = Vtmax_atAt, .solve_cstart_lstart = solve_cstart_lstart,
    .define_TRt = define_TRt
){
  # Variables
  c_path_normal = list()
  c_path_control = list()
  c_path_treat = list()
  l_path_normal = list()
  l_path_control = list()
  l_path_treat = list()
  A_path_normal = list()
  A_path_control = list()
  A_path_treat = list()
  
  
  find_controlvar_att = function (t, At, coeff, It, TRt, wt, lmax, σ, η, B, δ, qhyp, β) {
    ## find Anextt
    c(Vt, Anextt) %<-% .Vtmax_atAt(t, At, coeff, It, TRt, wt, lmax, σ, η, B, δ, qhyp, β)
    ## solve optimal bundle
    bundle = .solve_cstart_lstart(t, At, Anextt, It, TRt, wt, lmax, σ, η, B)
    cstart = bundle['cstart']
    lstart = bundle['lstart']
    return(c(lstart, cstart, Anextt))
  }
  
  # Initial asset   
  At_normal = A1
  At_control = A1
  At_treat = A1
  
  # Paths
  ## Normal
  for (t in seq(1, n_periods)) {
    horizon_period <- min((t + W - 1), length(coeffs_windows))
    coeffs <- coeffs_windows[[horizon_period]]
    wt <- w_bl*(1+.ω)^(t - .bl_period)
    TRt <- .define_TRt(t, 'normal', treat_transfer_monthly, control_transfer_monthly)
    c(lt_normal, ct_normal, Anextt_normal) %<-%  find_controlvar_att(t, At_normal, coeffs$normal, It, TRt, wt, lmax, σ, η, B, δ, qhyp, β)
    l_path_normal[t] <- lt_normal
    c_path_normal[t] <- ct_normal
    A_path_normal[t] <- At_normal
    At_normal <- Anextt_normal
  }
  
  
  ## Control
  for (t in seq(1, n_periods)) {
    horizon_period <- min((t + W - 1), length(coeffs_windows))
    coeffs <- coeffs_windows[[horizon_period]]
    wt <- w_bl*(1+.ω)^(t - .bl_period)
    TRt <- .define_TRt(t, 'control', treat_transfer_monthly, control_transfer_monthly)
    ### Should expect no change in TR before start of benefit
    if (t < .bi_start_period) {coeff_control_att = coeffs$normal} else {coeff_control_att = coeffs$control}
    c(lt_control, ct_control, Anextt_control) %<-%  find_controlvar_att(t, At_control, coeff_control_att, It, TRt, wt, lmax, σ, η, B, δ, qhyp, β)
    l_path_control[t] <- lt_control
    c_path_control[t] <- ct_control
    A_path_control[t] <- At_control
    At_control <- Anextt_control
  }
  
  ## Treat
  for (t in seq(1, n_periods)) {
    horizon_period <- min((t + W - 1), length(coeffs_windows))
    coeffs <- coeffs_windows[[horizon_period]]
    wt <- w_bl*(1+.ω)^(t - .bl_period)
    TRt <- .define_TRt(t, 'treat', treat_transfer_monthly, control_transfer_monthly)
    ### Should expect no change in TR before start of benefit
    if (t < .bi_start_period) {coeff_treat_att = coeffs$normal} else {coeff_treat_att = coeffs$treat}
    c(lt_treat, ct_treat, Anextt_treat) %<-%  find_controlvar_att(t, At_treat, coeff_treat_att, It, TRt, wt, lmax, σ, η, B, δ, qhyp, β)
    l_path_treat[t] <- lt_treat
    c_path_treat[t] <- ct_treat
    A_path_treat[t] <- At_treat
    At_treat <- Anextt_treat
  }
  
  return(list(c_path_normal, c_path_control, c_path_treat, l_path_normal, l_path_control, l_path_treat, A_path_normal, A_path_control, A_path_treat))
}


# 🎨 Plot ----
plot_paths = function (
    paths, start_period, end_period,
    .ω = ω,
    .n_periods_year = n_periods_year, .retirement_period = retirement_period,
    .bl_period = bl_period,
    .bi_start_period = bi_start_period, .bi_end_period = bi_end_period
) {
  paths = lapply(lapply(paths, '[', start_period:end_period), unlist)
  
  c_path_normal <- paths[[1]]
  c_path_control <- paths[[2]]
  c_path_treat <- paths[[3]]
  l_path_normal <- paths[[4]]
  l_path_control <- paths[[5]]
  l_path_treat <- paths[[6]]
  A_path_normal <- paths[[7]]
  A_path_control <- paths[[8]]
  A_path_treat <- paths[[9]]
  
  A_bl <- A_path_normal[[.bl_period]]
  
  wt = sapply(seq(1, end_period) - .bl_period, function(multiplier) w_bl*(1+.ω)^(multiplier))[start_period:end_period]
  
  
  
  fig <- ggplot() +
    geom_point(aes(x = seq(start_period:end_period), y = unlist(A_path_control) - A_bl, color = 'A_control'), alpha =.4) +
    geom_point(aes(x = seq(start_period:end_period), y = unlist(A_path_treat) - A_bl, color = 'A_treat'), alpha =.4) +
    geom_point(aes(x = seq(start_period:end_period), y = unlist(l_path_control)*wt , color = 'l_control'), alpha =.4) +
    geom_point(aes(x = seq(start_period:end_period), y = unlist(l_path_treat)*wt , color = 'l_treat'), alpha =.4) +
    geom_point(aes(x = seq(start_period:end_period), y = unlist(c_path_control) , color = 'c_control'), alpha =.4) +
    geom_point(aes(x = seq(start_period:end_period), y = unlist(c_path_treat), color = 'c_treat'), alpha =.4) +
    {if (start_period <= .bi_start_period) geom_vline(aes(xintercept = .bi_start_period - start_period + 1), color = 'cyan')} +
    {if (start_period <= .bi_end_period) geom_vline(aes(xintercept=.bi_end_period - start_period + 1), color='darkolivegreen3') }+
    {if (end_period >= .retirement_period) geom_vline(aes(xintercept=.retirement_period - start_period + 1), color='red')}+
    scale_x_continuous(labels = function(x) (x - .bi_start_period + start_period -1)/.n_periods_year)+ 
    scale_color_manual(values = c('A_control' = 'black', 'A_treat' = 'blue', 
                                  'l_control' = 'yellow', 'l_treat' = 'orange', 
                                  'c_control' = 'green', 'c_treat' = 'purple'),
                       labels = c('A_control' = 'Asset Control', 'A_treat' = 'Asset Treatment', 
                                  'l_control' = 'Labor Supply Control', 'l_treat' = 'Labor Supply Treatment', 
                                  'c_control' = 'Consumption Control', 'c_treat' = 'Consumption Treatment')) +
    theme_bw()+
    theme(
      axis.title = element_text(size = 20),  # Axis titles
      axis.text = element_text(size = 16),   # Axis text
      legend.text = element_text(size = 16),  # Legend text
      legend.title = element_blank()
    )+ 
    
    labs(x = "Years since start of program", y = "US Dollar")
  return(fig)
}


plot_effects <- function(
    paths, start_period, end_period, momentstype = "nominal", l_as_income = 0, paper = 0, admin = 0,
    .c_effect = c_effect, .l_effect = l_effect, .A_effect = A_effect,
    .c_effect_paper = c_effect_paper, .l_effect_paper = l_effect_paper,
    .income_effect_admin = income_effect_admin,
    .bl_period = bl_period, .ml_period = ml_period, .el_period = el_period,
    .bi_start_period = bi_start_period, .n_periods_year = n_periods_year,
    .ω = ω,
    .assemble_moments = assemble_moments, .assemble_weights = assemble_weights
){
  c_path_normal = unlist(paths[[1]])
  c_path_control = unlist(paths[[2]])
  c_path_treat = unlist(paths[[3]])
  l_path_normal = unlist(paths[[4]])
  l_path_control = unlist(paths[[5]])
  l_path_treat = unlist(paths[[6]])
  A_path_normal = unlist(paths[[7]])
  A_path_control = unlist(paths[[8]])
  A_path_treat = unlist(paths[[9]])
  wt = sapply(seq(1, length(l_path_treat)) - .bl_period, function(multiplier) w_bl*(1+.ω)^(multiplier))
  
  
  # Model Effect
  c_model_effect = (c_path_treat - c_path_control)
  l_model_effect = (l_path_treat - l_path_control)
  
  if (momentstype == "pct") {
    c_model_effect = c_model_effect/c_path_control
    l_model_effect = l_model_effect/l_path_control
  }

  
  # Effect from data
  if (momentstype == "pct") {
    .c_effect[-c(1,2)] =  .c_effect[-c(1,2)]/.c_effect$controlmean
    .l_effect[-c(1,2)] =  .l_effect[-c(1,2)]/.l_effect$controlmean
    .A_effect[-c(1,2)] =  .A_effect[-c(1,2)]/.A_effect$controlmean
    .income_effect_admin[-c(1,2)] =  .income_effect_admin[-c(1,2)]/.income_effect_admin$controlmean
  }
  
  effect_df = as_tibble(cbind(period=seq(1, length(c_path_normal)), c_model_effect = c_model_effect, l_model_effect = l_model_effect)) 
  
  c_real_effect <- .c_effect %>%
    mutate(period = round(.bi_start_period + month/12*.n_periods_year - 1)) %>% 
    select(period, beta, std.error, period) %>%
    rename(c_real_effect = beta,
           c_real_std.error = std.error)
  
  c_real_effect_paper <- .c_effect_paper %>%
    mutate(period = round(.bi_start_period + month/12*.n_periods_year - 1)) %>% 
    select(period, beta, std.error, period) %>%
    rename(c_real_effect_paper = beta,
           c_real_std.error_paper = std.error)
  
  
  l_real_effect <- .l_effect %>%
    mutate(period = round(.bi_start_period + month/12*.n_periods_year - 1)) %>%
    select(period, beta, std.error, period) %>%
    rename(l_real_effect = beta,
           l_real_std.error = std.error)
  
  l_real_effect_paper <- .l_effect_paper %>%
    mutate(period = round(.bi_start_period + month/12*.n_periods_year - 1)) %>%
    select(period, beta, std.error, period) %>%
    rename(l_real_effect_paper = beta,
           l_real_std.error_paper = std.error)
  
  income_real_effect <- .income_effect_admin %>%
    mutate(period = round(.bi_start_period + month/12*.n_periods_year - 1)) %>%
    select(period, beta, std.error, period) %>%
    rename(income_real_effect = beta,
           income_real_std.error = std.error)
  
  
  effect_df = effect_df %>%
    merge(c_real_effect, all = T, by = "period") %>%
    merge(l_real_effect, all = T, by = "period") %>%
    merge(c_real_effect_paper, all = T, by = "period") %>%
    merge(l_real_effect_paper, all = T, by = "period") %>%
    merge(income_real_effect, all = T, by = "period") %>%
    mutate(wt = c(wt, c(0))) %>%
    {if (momentstype == "nominal") {
      mutate(.,
             l_real_effect = l_real_effect * wt,
             l_real_std.error = l_real_std.error*wt,
             l_real_effect_paper = l_real_effect_paper * wt,
             l_real_std.error_paper = l_real_std.error_paper*wt,
             l_model_effect = l_model_effect * wt
             )
    } else .}
  
  effect_df %>%
    filter(period >= start_period & period <= end_period) %>%
    ggplot(mapping = aes(x = period)) +
    geom_line(aes(y = c_model_effect, color = "E1"), alpha = 0.5) +
    geom_point(aes(y = c_model_effect, color = "E1"), alpha = 0.5) +
    {if(paper == 0) geom_point(aes(y = c_real_effect, color = "E2"))}+
    {if(paper == 1) geom_point(aes(y = c_real_effect_paper, color = "E2"))}+
    geom_line(aes(y = l_model_effect, color = "E3"), alpha = 0.5) +
    geom_point(aes(y = l_model_effect, color = "E3"), alpha = 0.5) +
    {if(paper == 0 & admin == 0) geom_point(aes(y = l_real_effect, color = "E4"))} +
    {if(paper == 1 & admin == 0) geom_point(aes(y = l_real_effect_paper, color = "E4"))} +
    {if(admin == 1) geom_point(aes(y = income_real_effect, color = "E4"))} +
    scale_x_continuous(labels = function(x) (x - .bi_start_period + 1)/.n_periods_year) +
    scale_color_manual(
      values = c("E1" = "green", "E2" = "blue", "E3" = "orange", "E4" = "red"),
      labels = c("E1" = "Consumption Model ", "E2" = "Consumption Observed", "E3" = "Labor Supply Model", "E4" = "Labor Supply Observed")
    ) + 
    theme_bw() +
    theme(
      axis.title = element_text(size = 20),  # Axis titles
      axis.text = element_text(size = 16),   # Axis text
      legend.text = element_text(size = 16),  # Legend text
      legend.title = element_blank()
    )+ 

    labs(y = "Treatment Effect", x = "Years since start of program")
}

# 🌸 CALCULATE MPE ----
calculate_MPE = function(
    paths, It, dI,
    .bl_period = bl_period, .ml_period = ml_period, .el_period = el_period, .ω = ω,
    .bi_start_period = bi_start_period, .bi_end_period = bi_end_period
){
  c_path_normal = unlist(paths[[1]])
  c_path_control = unlist(paths[[2]])
  c_path_treat = unlist(paths[[3]])
  l_path_normal = unlist(paths[[4]])
  l_path_control = unlist(paths[[5]])
  l_path_treat = unlist(paths[[6]])
  A_path_normal = unlist(paths[[7]])
  A_path_control = unlist(paths[[8]])
  A_path_treat = unlist(paths[[9]])
  wt = sapply(seq(1, length(l_path_treat)) - .bl_period, function(multiplier) w_bl*(1+.ω)^(multiplier))
  
  
  
  dI_temp = rep(0, length(A_path_control))
  dI_temp[.bi_start_period: .bi_end_period] <- rep(dI, (.bi_end_period - .bi_start_period + 1))
  dI <- dI_temp
  
  
  l_model_effect = (l_path_treat - l_path_control)
  inc_hh_effect = l_model_effect*wt
  savings_effect = (A_path_treat) - (A_path_control)
  delta_savings_effect = savings_effect - lag(savings_effect)

  mpe = inc_hh_effect/(dI - delta_savings_effect)
  
  mpe_ml = mpe[.ml_period]
  mpe_el = mpe[.el_period]
  mpe_pool = mpe_ml*0.3 + mpe_el*0.7
  
  
  return(c(mpe_ml = mpe_ml, mpe_el = mpe_el, mpe_pool = mpe_pool))
  
}


# 🏋WEIGHTED DISTANCE ----
assemble_moments = function(
    momentsversion, momentstype,
    .c_effect = c_effect, .l_effect = l_effect, .A_effect = A_effect,
    .c_effect_paper = c_effect_paper, .l_effect_paper = l_effect_paper,
    .income_effect_admin = income_effect_admin
) {
  
  if (momentstype == "pct") {
    .c_effect[-c(1,2)] =  .c_effect[-c(1,2)]/.c_effect$controlmean
    .c_effect_paper[-c(1,2)] =  .c_effect_paper[-c(1,2)]/.c_effect_paper$controlmean
    .l_effect[-c(1,2)] =  .l_effect[-c(1,2)]/.l_effect$controlmean
    .l_effect_paper[-c(1,2)] = .l_effect_paper[-c(1,2)]/.l_effect_paper$controlmean
    .A_effect[-c(1,2)] =  .A_effect[-c(1,2)]/.A_effect$controlmean
    .income_effect_admin[-c(1,2)] =  .income_effect_admin[-c(1,2)]/.income_effect_admin$controlmean
  }
  
  if (momentsversion == 1) data_moments  = .l_effect$beta
  if (momentsversion == 2) data_moments  = c(.c_effect$beta, .l_effect$beta)
  if (momentsversion == 3) data_moments  = c(.c_effect$beta, .l_effect$beta, .A_effect$value)
  if(momentsversion == 4) data_moments  = c(.c_effect_paper$beta, .l_effect_paper$beta)
  if(momentsversion == 5) data_moments  = c(.c_effect_paper$beta, .l_effect_paper$beta, .A_effect$value)
  if(momentsversion == 6) data_moments  = c(.c_effect$beta, .income_effect_admin$beta, .A_effect$value)
  if(momentsversion == 7) data_moments  = c(.c_effect$beta, .income_effect_admin$beta)
  
  return(data_moments)
  
}


assemble_weights = function(
    momentsversion, momentstype,
    .c_effect = c_effect, .l_effect = l_effect, 
    .c_effect_paper = c_effect_paper, .l_effect_paper = l_effect_paper,
    .A_effect = A_effect, .income_effect_admin = income_effect_admin
){
  if (momentstype == "pct") {
    .c_effect[-c(1,2)] =  .c_effect[-c(1,2)]/.c_effect$controlmean
    .c_effect_paper[-c(1,2)] =  .c_effect_paper[-c(1,2)]/.c_effect_paper$controlmean
    .l_effect[-c(1,2)] =  .l_effect[-c(1,2)]/.l_effect$controlmean
    .l_effect_paper[-c(1,2)] = .l_effect_paper[-c(1,2)]/.l_effect_paper$controlmean
    .A_effect[-c(1,2)] =  .A_effect[-c(1,2)]/.A_effect$controlmean
    .income_effect_admin[-c(1,2)] =  .income_effect_admin[-c(1,2)]/.income_effect_admin$controlmean
  }
  
  c_weights = (.c_effect$std.error * ((1/3*1/nrow(.c_effect))^(-1/2)))**2
  c_weights_paper = (.c_effect_paper$std.error *1/3*c(0.2, 0.3, 0.5))**2
  l_weights = (.l_effect$std.error * ((1/3*1/nrow(.l_effect))^(-1/2)))**2
  l_weights_paper = (.l_effect_paper$std.error *1/3*c(0.3, 0.7))**2
  A_weights = (.A_effect$std.error * ((1/3*1/nrow(.A_effect))^(-1/2)))**2
  A_weights_paper = (.A_effect$std.error * ((1/3*1/2*c(0.3, 0.7, 0.3, 0.7))^(-1/2)))**2
  income_weights = (.income_effect_admin$std.error * ((1/3*1/nrow(.income_effect_admin))^(-1/2)))**2
  
  if (momentsversion == 1) weight_matrix  = diag(l_weights)
  if (momentsversion == 2) weight_matrix  = diag(c(c_weights, l_weights))
  if (momentsversion == 3) weight_matrix  = diag(c(c_weights, l_weights, A_weights))
  if (momentsversion == 4) weight_matrix = diag(c(c_weights_paper, l_weights_paper))
  if (momentsversion == 5) weight_matrix = diag(c(c_weights_paper, l_weights_paper, A_weights_paper))
  if (momentsversion == 6) weight_matrix = diag(c(c_weights, income_weights, A_weights))
  if (momentsversion == 7) weight_matrix = diag(c(c_weights, income_weights))
  
  
  return(weight_matrix)
}


calculate_weighted_distance = function(
    paths, momentsversion = 1, momentstype = "pct",
    .ω = ω, .w_bl = w_bl,
    .c_effect = c_effect, .l_effect = l_effect, .A_effect = A_effect,
    .c_effect_paper = c_effect_paper, .l_effect_paper = l_effect_paper,
    .income_effect_admin = income_effect_admin,
    .bl_period = bl_period, .ml_period = ml_period, .el_period = el_period,
    .bi_start_period = bi_start_period, .n_periods_year = n_periods_year,
    .assemble_moments = assemble_moments, .assemble_weights = assemble_weights
){
  
  data_moments <- .assemble_moments(momentsversion, momentstype)
  weight_matrix <- .assemble_weights(momentsversion, momentstype)
  
  
  
  c_path_normal = unlist(paths[[1]])
  c_path_control = unlist(paths[[2]])
  c_path_treat = unlist(paths[[3]])
  l_path_normal = unlist(paths[[4]])
  l_path_control = unlist(paths[[5]])
  l_path_treat = unlist(paths[[6]])
  A_path_normal = unlist(paths[[7]])
  A_path_control = unlist(paths[[8]])
  A_path_treat = unlist(paths[[9]])
  wt = sapply(seq(1, length(l_path_treat)) - .bl_period, function(multiplier) .w_bl*(1+.ω)^(multiplier))
  income_path_control = l_path_control * wt
  
  
  c_model_effect = (c_path_treat - c_path_control)
  l_model_effect = (l_path_treat - l_path_control)
  A_model_effect = (A_path_treat - A_path_control)
  income_model_effect = l_model_effect * wt
  
  if (momentstype == "pct") {
    c_model_effect = c_model_effect/c_path_control
    l_model_effect = l_model_effect/c_path_control
    income_model_effect = income_model_effect/income_path_control
  }
  
  c_effect <- .c_effect %>%
    mutate(period = round(.bi_start_period + month/12*.n_periods_year - 1))
  
  c_effect_paper <- .c_effect_paper %>%
    mutate(period = round(.bi_start_period + month/12*.n_periods_year - 1))
  
  l_effect <- .l_effect %>%
    mutate(period = round(.bi_start_period + month/12*.n_periods_year - 1))
  
  l_effect_paper <- .l_effect_paper %>%
    mutate(period = round(.bi_start_period + month/12*.n_periods_year - 1))
  
  A_effect <- .A_effect %>%
    mutate(period = round(.bi_start_period + month/12*.n_periods_year - 1))
  
  income_effect_admin <- .income_effect_admin %>%
    mutate(period = round(.bi_start_period + month/12*.n_periods_year - 1))
  
  c_model_moments <- c_model_effect[c_effect$period]
  c_model_moments_paper <- c_model_effect[c_effect_paper$period]
  l_model_moments <- l_model_effect[l_effect$period]
  l_model_moments_paper <- l_model_effect[l_effect_paper$period]
  A_model_moments_control_raw = A_path_control[A_effect$period[1:2]]
  A_model_moments_treat_effect = A_model_effect[A_effect$period[3:4]]
  A_model_moments <- c(A_model_moments_control_raw, A_model_moments_treat_effect)
  income_model_moments_admin <- income_model_effect[income_effect_admin$period]
  
  if (momentsversion == 1) model_moments <- l_model_moments
  if (momentsversion == 2) model_moments <- c(c_model_moments, l_model_moments)
  if (momentsversion == 3) model_moments <- c(c_model_moments, l_model_moments, A_model_moments)
  if (momentsversion == 4) model_moments <- c(c_model_moments_paper, l_model_moments_paper)
  if (momentsversion == 5) model_moments <- c(c_model_moments_paper, l_model_moments_paper, A_model_moments)
  if (momentsversion == 6) model_moments <- c(c_model_moments, income_model_moments_admin, A_model_moments)
  if (momentsversion == 7) model_moments <- c(c_model_moments, income_model_moments_admin)
  
  
  distance = t(model_moments - data_moments) %*% pracma::inv(weight_matrix) %*%  (model_moments - data_moments)
  
  if ((momentsversion == 3 | momentsversion == 4| momentsversion == 5| momentsversion == 6 | momentsversion == 7) & momentstype == "pct") distance = c(NA)
  
  return(distance[[1,1]])
}
