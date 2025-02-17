r_yearly = 0.02
r_borrowing_yearly = 0.1
c_bl_yearly = c(3312.499, 1607.245) * 12 # control mean and sd
c_ml_yearly = c(4222.650, 1957.721) * 12
c_bl_yearly = c_bl_yearly
c_ml_yearly = c_ml_yearly
median_hh_income_yearly = 76660 # from St Louis FRED, personal income
γc_pct <- 0.12
γc_yearly <-median_hh_income_yearly * γc_pct
l_bl_yearly = c(30.2, 26.5) *52 # joint_work_hrs
l_ml_yearly = c(39.1, 27.9) * 52
w_bl_pretax = c(13.0, 6.61)
w_ml_pretax = c(16.8, 9.25)
mean_hh_ATR <- 0.131
ω_yearly = 0.06 # %Δwage between el & ml - real interest rate 
w_bl <- w_ml_pretax[1]*(1+ω_yearly)**(-1.5) * (1-mean_hh_ATR) # ml is more accurate
w_ml <- w_ml_pretax[1]* (1-mean_hh_ATR)
A_bl = -4672.422 
It_yearly = 43325.72*0.02 # control mean sd hh_inc_total * 0.02 
treat_transfer_monthly <- 1000
control_transfer_monthly <- 50
lmax_yearly <- 60*52*2



r <- (1+r_yearly)^(1/n_periods_year)-1
r_borrowing <- (1+r_borrowing_yearly)^(1/n_periods_year)-1
c_bl <- c_bl_yearly/n_periods_year
c_ml <- c_ml_yearly/n_periods_year
γc <- γc_yearly/n_periods_year
l_bl <- l_bl_yearly/n_periods_year
l_ml <- l_ml_yearly/n_periods_year
lmax <- lmax_yearly/n_periods_year
It = It_yearly/n_periods_year
ω <- (1+ω_yearly)^(1/n_periods_year)-1
A_min = min(A_bl, 0) # Assume that agent have borrow the maximum they can borrow at BL










