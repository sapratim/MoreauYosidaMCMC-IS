# Run from the TrendFiltering directory:
source("Timings_check/timings_functions.R")
load("warmup.Rdata")

iter <- 2000
L <- 100
lamb_coeff <- 0.001
alpha_hat <- 5
sigma2_hat <- 9
eps_my <- 0.015
eps_px <- 0.00035
tol <- 1e-5
max_iter <- 200L
trend_control <- trendfilter.control.list(obj_tol = tol, max_iter = max_iter)
D_mat <- getD(k = 1, n = length(y), x = x)


set.seed(54321)
my <- myhmc_record(y, alpha_hat, sigma2_hat, k = 1, grid = x,
                   iter = iter, eps_hmc = eps_my, L = L, start = warmup_end_iter)
set.seed(54321)
px <- pxhmc_record(y, alpha_hat, sigma2_hat, k = 1, grid = x,
                   iter = iter, eps_hmc = eps_px, L = L, start = warmup_end_iter)

# Saves both final-proposal records and all intermediate leapfrog ADMM/roughness records.
save(iter, my, px, file = "Timings_check/small_run.Rdata")
