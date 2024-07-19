
source("TF_functions.R")
load("warmup.Rdata")

iter <- 1e5
lamb_coeff <- 0.001
delta_samp_is <- 0.0015
delta_samp_pxm <- 0.0008
delta_bark_my <- 0.0014
delta_bark_px <- 0.0009
D_mat <- getD(k=1, n=1e2, x)   #  D matrix
alpha_hat <- 5   # obtained from the first dataset
sigma2_hat <- 9  # obtained from the first dataset
k <- 1

mala.is <- mymala(y, alpha_hat, sigma2_hat, k=1, grid=x, iter = iter, 
                  delta = delta_samp_is, start = warmup_end_iter)

pxmala.run <- px.mala(y, alpha_hat, sigma2_hat, k=1, grid=x, iter = iter, 
                      delta = delta_samp_pxm, start = warmup_end_iter)

mybark <- mybarker(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter,
                   delta = delta_bark_my, start = warmup_end_iter)

pxbark <- px.barker(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter,
                   delta = delta_bark_px, start = warmup_end_iter)

my.hmc <- myhmc(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter,
                eps_hmc = 0.015, L = 100, start = warmup_end_iter)

px.hmc <- pxhmc(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter,
                eps_hmc = 0.0003, L = 100, start = warmup_end_iter) 

output_single_run_mala <- list(mala.is[[1]], pxmala.run)
output_single_run_bark <- list(mybark[[1]],  pxbark[[1]])
output_single_run_hmc <- list(my.hmc[[1]], px.hmc[[1]])
log_wts <- list(mala.is[[2]], mybark[[2]], my.hmc[[2]])

save(output_single_run_mala, file = "single_chain_mala.Rdata")
save(output_single_run_bark, file = "single_chain_bark.Rdata")
save(output_single_run_hmc, file = "single_chain_hmc.Rdata")
save(log_wts, file = "single_chain_log_weights.Rdata")
