# 🌝 SETUP ----
setwd("/Users/leodai/Library/CloudStorage/GoogleDrive-lckdai716@gmail.com/My Drive/Eva/Basic Income/dynamic_labor_supply")
# install.packages("pacman")
pacman::p_load("dplyr", "tidyr", "stringr", "furrr", "purrr", "future", "gtools", 'zeallot', 'DEoptim', 'plotly', 'datawizard', 'progressr', 'expm', 'npreg', 'pracma', 'xtable')
plan(multisession, workers = detectCores() -1)
set.seed(167898)

# args = commandArgs(trailingOnly=TRUE)
# hhcomposition = as.integer(args[1])
hhcomposition = 4
# 1: single worker 2: dual workers, exogenous partner income
# 3: dual workers, unitary household # 4: weighted avg of 1 and 3
hhcomposition_dir <- list.files(".", sprintf("hhcomposition%s.*", hhcomposition),include.dirs = T)

# momentsversion
# 1: labor supply + consumption 2: statusquo 3: labor supply + consumption + asset 

# momentstypd
# 1: pct 2: nominal




# 📦 MODEL FUNCTIONS ----
source("dynamic_labor_supply_functions.R")

# 🌝 MODEL PARAMETERS ----
source("dynamic_labor_supply_model_parameters.R")


# 🌝 ECONOMIC PARAMETERS ----
source(sprintf("%s/dynamic_labor_supply_economic_parameters.R", hhcomposition_dir))


# 💃 MOMENTS AND DISTANCE ----
source(sprintf("%s/dynamic_labor_supply_moments.R", hhcomposition_dir))





# 🌵 STANDARD MODEL ----
## optimization functions ----
find_σηβδ_distance <- function(σ, η, β, δ_yearly){
  print(sprintf("σ: %s, η: %s, β: %s, δ: %s", σ, η, β, δ_yearly))
  δ <- δ_yearly^(1/n_periods_year)
  B_grid = get_B_grid(σ, η, 1, c_ml, l_ml)
  
  
  distance_df <- tibble(
    "σ" = numeric(),
    "η" = numeric(),
    "β" = numeric(),
    "δ_yearly" = numeric(),
    "B" = numeric(),
    "distance1" = numeric(),
    "distance2" = numeric(),
    "distance3" = numeric(),
    "distance4" = numeric(),
    "distance5" = numeric(),
    "distance6" = numeric(),
    "distance7" = numeric(),
    "distance8" = numeric(),
    "distance9" = numeric()
    )
  
  for (B in B_grid) {
    models_bystatus = solve_dynamic_labor_supply_model(terminal_period, It, w_bl, lmax, σ, η, B, δ, treat_transfer_monthly, control_transfer_monthly)
    paths = get_paths(A_bl, bi_end_period + 2*n_periods_year, models_bystatus, w_bl, It, lmax, σ, η, B, δ, treat_transfer_monthly, control_transfer_monthly, qhyp = 1, β)
    
    distance1 = calculate_weighted_distance(paths, momentsversion = 1, momentstype = "pct")
    distance2 = calculate_weighted_distance(paths, momentsversion = 2, momentstype = "pct")
    distance3 = calculate_weighted_distance(paths, momentsversion = 1, momentstype = "nominal")
    distance4 = calculate_weighted_distance(paths, momentsversion = 2, momentstype = "nominal")
    distance5 = calculate_weighted_distance(paths, momentsversion = 3, momentstype = "nominal")
    distance6 = calculate_weighted_distance(paths, momentsversion = 4, momentstype = "nominal")
    distance7 = calculate_weighted_distance(paths, momentsversion = 5, momentstype = "nominal")
    distance8 = calculate_weighted_distance(paths, momentsversion = 6, momentstype = "nominal")
    distance9 = calculate_weighted_distance(paths, momentsversion = 6, momentstype = "nominal")
    print(sprintf("σ: %s, η: %s, β: %s, δ: %s, B: %s - %s, %s, %s, %s, %s, %s, %s, %s", σ, η, β, δ_yearly, B, distance1, distance2, distance3, distance4, distance5, distance6, distance7, distance8, distance9))
    
    distance_df <- distance_df %>%
      bind_rows(
        tibble(
          "σ" = σ,
          "η" = η,
          "β" = β,
          "δ_yearly" = δ_yearly,
          "B" = B,
          "distance1" = distance1,
          "distance2" = distance2,
          "distance3" = distance3,
          "distance4" = distance4,
          "distance5" = distance5,
          "distance6" = distance6,
          "distance7" = distance7,
          "distance8" = distance8,
          "distance9" = distance9
        )
      )
    
  }
    
  return(distance_df)
}


find_σηβδ_distance_wrapper <- function(params){
  σ = params[[1]]
  η = params[[2]]
  β = params[[3]]
  δ_yearly = params[[4]]
  
  distance_df <- find_σηβδ_distance(σ, η, β, δ_yearly)
  return(distance_df)
  
}

# ## optimize ----
σ_grid = c(1, 2)
η_grid = seq(0.05, 0.3, length.out = 6)
β_grid = seq(0.65, 1, length.out = 8)
δ_grid = seq(0.95, 1.05, length.out = 11)


σηβδ_grid = expand.grid(σ = σ_grid, η = η_grid, β = β_grid, δ= δ_grid)
σηβδ_distances = future_pmap(σηβδ_grid, find_σηβδ_distance)
σηβδ_distance <- do.call("rbind", σηβδ_distances)
saveRDS(σηβδ_distance, sprintf("%s/dynamic_labor_supply_standard_σηβδ_distance_variousmoments_Dec29.rds", hhcomposition_dir))


## optimized output ----
σηβδ_distance <- readRDS(sprintf("%s/dynamic_labor_supply_standard_σηβδ_distance_variousmoments_Dec29.rds", hhcomposition_dir))

r_yearly = 0.02
r <- (1+r_yearly)^(1/n_periods_year)-1



### Version 3 Nominal Best Fit ----
momentsversion = 3
momentstype = "nominal"

σηβδ_distance <- σηβδ_distance %>% arrange(distance5)
σ = σηβδ_distance %>% slice(1) %>% pull(σ)
η = σηβδ_distance %>% slice(1) %>% pull(η)
β = σηβδ_distance %>% slice(1) %>% pull(β)
δ_yearly = σηβδ_distance %>% slice(1) %>% pull(δ_yearly)

δ <- δ_yearly^(1/n_periods_year)
B = compute_B(σ, η, c_ml[[1]], l_ml[[1]], w_ml, γc)

models_bystatus = solve_dynamic_labor_supply_model(terminal_period, It, w_bl, lmax, σ, η, B, δ, treat_transfer_monthly, control_transfer_monthly)
paths = get_paths(A_bl, terminal_period, models_bystatus, w_bl, It, lmax, σ, η, B, δ, treat_transfer_monthly, control_transfer_monthly, qhyp = 1, β)

fig_path2 <- plot_paths(paths, 1, bi_end_period + 2*n_periods_year)
fig_moments2nominal <- plot_effects(paths, bi_start_period -1 , bi_end_period + 2*n_periods_year, "nominal")

ggsave(sprintf("%s/figure_paths_biplus2yrs_momentscompv%s_momentstype%s_bestfit.pdf", hhcomposition_dir, momentsversion, momentstype), plot = fig_path2, width = 10, height = 10)
ggsave(sprintf("%s/figure_moments_biplus2yrs_momentscompv%s_momentstype%s_against%s_bestfit.pdf", hhcomposition_dir, momentsversion, momentstype, "nominal"), plot = fig_moments2nominal, width = 10, height = 10)


MPE5 <- calculate_MPE(paths, It, dI = 950*12/n_periods_year)
params5 <- c("σ" = σ, "η" = η, "β" = β, "δ_yearly" = δ_yearly, "r_yearly" = r_yearly)



### Version 3 Nominal rate Golosov----
momentsversion = 3
momentstype = "nominal"

r_yearly = 0.025
δ_yearly = 1/(1+0.025)

σηβδ_distance <- σηβδ_distance %>% arrange(distance5)
σ = σηβδ_distance %>% slice(1) %>% pull(σ)
η = σηβδ_distance %>% slice(1) %>% pull(η)
β = σηβδ_distance %>% slice(1) %>% pull(β)

r <- (1+r_yearly)^(1/n_periods_year)-1
δ <- δ_yearly^(1/n_periods_year)
B = compute_B(σ, η, c_ml[[1]], l_ml[[1]], w_ml, γc)

models_bystatus = solve_dynamic_labor_supply_model(terminal_period, It, w_bl, lmax, σ, η, B, δ, treat_transfer_monthly, control_transfer_monthly)
paths = get_paths(A_bl, terminal_period, models_bystatus, w_bl, It, lmax, σ, η, B, δ, treat_transfer_monthly, control_transfer_monthly, qhyp = 1, β)



fig_path2 <- plot_paths(paths, 1, bi_end_period + 2*n_periods_year)
fig_moments2nominal <- plot_effects(paths, bi_start_period -1 , bi_end_period + 2*n_periods_year, "nominal")
ggsave(sprintf("%s/figure_paths_biplus2yrs_momentscompv%s_momentstype%s_rateGolosov.pdf", hhcomposition_dir, momentsversion, momentstype), plot = fig_path2, width = 10, height = 10)
ggsave(sprintf("%s/figure_moments_biplus2yrs_momentscompv%s_momentstype%s_against%s_rateGolosov.pdf", hhcomposition_dir, momentsversion, momentstype, "nominal"), plot = fig_moments2nominal, width = 10, height = 10)

MPE53 <- calculate_MPE(paths, It, dI = 950*12/n_periods_year)
params53 <- c("σ" = σ, "η" = η, "β" = β, "δ_yearly" = δ_yearly, "r_yearly" = r_yearly)



### Version 6 with Admin data ----
momentsversion = 6
momentstype = "nominal"

r_yearly = 0.025
δ_yearly = 1/(1+0.025)

σηβδ_distance <- σηβδ_distance %>% arrange(distance8) %>% filter(σ != 1.5)
σ = σηβδ_distance %>% slice(1) %>% pull(σ)
η = σηβδ_distance %>% slice(1) %>% pull(η)
β = σηβδ_distance %>% slice(1) %>% pull(β)
δ_yearly = σηβδ_distance %>% slice(1) %>% pull(δ_yearly)

δ <- δ_yearly^(1/n_periods_year)
B = compute_B(σ, η, c_ml[[1]], l_ml[[1]], w_ml, γc)

models_bystatus = solve_dynamic_labor_supply_model(terminal_period, It, w_bl, lmax, σ, η, B, δ, treat_transfer_monthly, control_transfer_monthly)
paths = get_paths(A_bl, terminal_period, models_bystatus, w_bl, It, lmax, σ, η, B, δ, treat_transfer_monthly, control_transfer_monthly, qhyp = 1, β)


fig_path2 <- plot_paths(paths, 1, bi_end_period + 2*n_periods_year)
fig_moments2nominal <- plot_effects(paths, bi_start_period -1 , bi_end_period + 2*n_periods_year, "nominal", admin = 1)

ggsave(sprintf("%s/figure_paths_biplus2yrs_momentscompv%s_momentstype%s_withadmin_bestfit.pdf", hhcomposition_dir, momentsversion, momentstype), plot = fig_path2, width = 10, height = 10)
ggsave(sprintf("%s/figure_moments_biplus2yrs_momentscompv%s_momentstype%s_withadmin_against%s_bestfit.pdf", hhcomposition_dir, momentsversion, momentstype, "nominal"), plot = fig_moments2nominal, width = 10, height = 10)


MPE8 <- calculate_MPE(paths, It, dI = 950*12/n_periods_year)
params8 <- c("σ" = σ, "η" = η, "β" = β, "δ_yearly" = δ_yearly, "r_yearly" = r_yearly)


### Version 6 Nominal rate Golosov----
momentsversion = 6
momentstype = "nominal"

r_yearly = 0.025
δ_yearly = 1/(1+0.025)

σηβδ_distance <- σηβδ_distance %>% arrange(distance8) %>% filter(σ != 1.5)
σ = σηβδ_distance %>% slice(1) %>% pull(σ)
η = σηβδ_distance %>% slice(1) %>% pull(η)
β = σηβδ_distance %>% slice(1) %>% pull(β)

r <- (1+r_yearly)^(1/n_periods_year)-1
δ <- δ_yearly^(1/n_periods_year)
B = compute_B(σ, η, c_ml[[1]], l_ml[[1]], w_ml, γc)

models_bystatus = solve_dynamic_labor_supply_model(terminal_period, It, w_bl, lmax, σ, η, B, δ, treat_transfer_monthly, control_transfer_monthly)
paths = get_paths(A_bl, terminal_period, models_bystatus, w_bl, It, lmax, σ, η, B, δ, treat_transfer_monthly, control_transfer_monthly, qhyp = 1, β)



fig_path2 <- plot_paths(paths, 1, bi_end_period + 2*n_periods_year)
fig_moments2nominal <- plot_effects(paths, bi_start_period -1 , bi_end_period + 2*n_periods_year, "nominal", admin = 1)
ggsave(sprintf("%s/figure_paths_biplus2yrs_momentscompv%s_momentstype%s_withadmin_rateGolosov.pdf", hhcomposition_dir, momentsversion, momentstype), plot = fig_path2, width = 10, height = 10)
ggsave(sprintf("%s/figure_moments_biplus2yrs_momentscompv%s_momentstype%s_withadmin_against%s_rateGolosov.pdf", hhcomposition_dir, momentsversion, momentstype, "nominal"), plot = fig_moments2nominal, width = 10, height = 10)

MPE83 <- calculate_MPE(paths, It, dI = 950*12/n_periods_year)
params83 <- c("σ" = σ, "η" = η, "β" = β, "δ_yearly" = δ_yearly, "r_yearly" = r_yearly)


### Output ----

mpe_table <- tribble(
  ~"Model", ~"\\sigma", ~"\\eta", ~"\\beta", ~"\\delta^{\\text{yearly}}", ~"r^{\\text{yearly}}", ~"MPE^{\\text{ml}}", ~"MPE^{\\text{el}}", ~"MPE^{\\text{pool}}",
  "Observed estimates", NA, NA, NA, NA, NA, eval(-2893.6*(-623.0982/-2070.2)/(950*12)), eval(-4742.2*(-1430/-2452.2)/(950*12)), eval(-4124.7*(-1503.1/-2504.3)/(950*12)),
  "Model - all estimates", params8[1], params8[2], params8[3], params8[4], params8[5], MPE8[1], MPE8[2],MPE8[3],
  "Model - Golosov et. al r and \\delta", params83[1], params83[2], params83[3], params83[4], params83[5], MPE83[1], MPE83[2], MPE83[3], 
)

saveRDS(mpe_table, sprintf("%s/mpe_table.rds", hhcomposition_dir))

mpe_table_tex <- mpe_table %>%
  
  mutate(across(colnames(.)[2:4], ~round(., 2))) %>%
  mutate(across(colnames(.)[5], ~round(., 3))) %>%
  mutate(across(matches("MPE"), ~round(., 2))) %>%
  mutate_all(as.character) %>%
  
  add_row(!!!setNames(as.list(colnames(.)), colnames(.)), .before = 1) %>%
  mutate_all(~replace_na(.x, "NA")) %>%
  mutate_all( ~paste0(.x, " & ")) %>%
  add_row(!!!setNames(rep("", length(colnames(.))), colnames(.)), .before=2) %>%
  mutate(Model = ifelse(row_number() == 2, "\\hline", Model)) %>%
  mutate(next_line = "\\\\") 

mpe_table_tex

mpe_table_tex %>%
  write.table(, file = sprintf("%s/mpe_table.tex", hhcomposition_dir), sep = "", row.names = FALSE, col.names = FALSE, quote = FALSE)

mpe_table



















