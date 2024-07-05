source("IS_trendf_functions_Pereyra.R")
load("covmat.Rdata")
load("MC_pcm.Rdata")

iter <- 1e5
delta_samp_is <- 0.095
delta_samp_pxm <- 0.015
lamb_coeff <- 0.001
D_mat <- getD(k=1, n=1e2, x)   #  D matrix

mala.is <- mymala(y, alpha_hat, sigma2_hat, k=1, grid=x, iter = iter, 
                  delta = delta_samp_is, covmat = covmat)

pxmala.run <- px.mala(y, alpha_hat, sigma2_hat, k=1, grid=x, iter = iter, 
                      delta = delta_samp_pxm, covmat = covmat)

my.hmc <- myhmc(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter,
                eps_hmc = 0.015, L = 100)

px.hmc <- pxhmc(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter,
                eps_hmc = 0.0003, L = 100) 
output_single_run <- list(mala.is[[1]], pxmala.run, my.hmc[[1]], px.hmc[[1]])
save(output_single_run, file = "single_chain.Rdata")
