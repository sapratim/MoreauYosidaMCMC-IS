rm(list = ls())
set.seed(8024248)
source("nuclear_norm_functions.R")
load("start_value.Rdata")
iter <- 1e5
lamb_coeff <- 1e-4
sigma2_hat <- 0.01
alpha_hat <- 1.15/sigma2_hat
step_ismala <- 0.00012
step_pxmala <- 0.0001
eps_is <- 0.008
eps_px <-  0.004
L <- 10

result_pxm <- px.mala(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat,
                      iter = iter, delta = step_pxmala, start = start_value)

result_is <- mymala(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat,
                    iter = iter, delta = step_ismala, start = start_value)

result_pxhmc <- pxhmc(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat, 
                      iter = iter, eps_hmc = eps_px, L=L, start = start_value)

result_ishmc <- myhmc(y=y, alpha = alpha_hat,lambda = lamb_coeff, sigma2 = sigma2_hat, 
                      iter = iter, eps_hmc = eps_is, L=L, start = start_value)

output_single <- list(result_is[[1]], result_pxm, result_ishmc[[1]], result_pxhmc[[1]])

save(output_single, file = "output_single_chain.Rdata")
