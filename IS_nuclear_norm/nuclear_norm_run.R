rm(list = ls())
source("nuclear_norm_functions.R")
load("start_value.Rdata")
iter <- 1e5 
lamb_coeff <- 0.00001
sigma2_hat <- 0.01
alpha_hat <- 1.15/sigma2_hat
step_ismala <- 0.00008
step_pxmala <- 0.00008

result_is <- mymala(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat,
                      iter = iter, delta = step_ismala, start = start_value)
result_pxm <- px.mala(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat,
                      iter = iter, delta = step_pxmala, start = start_value)
result_ishmc <- myhmc(y=y, alpha = alpha_hat,lambda = lamb_coeff, sigma2 = sigma2_hat, 
                      iter = iter, eps_hmc = 0.0047, L=10, start = start_value)
result_pxhmc <- pxhmc(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat, 
                      iter = iter, eps_hmc = 0.0047, L=10, start = start_value)

save(result_is, file = "IS_MALA.Rdata")
save(result_pxm, file = "Px_MALA.Rdata")
save(result_ishmc, file = "IS_HMC.Rdata")
save(result_pxhmc, file = "Px_HMC.Rdata")
