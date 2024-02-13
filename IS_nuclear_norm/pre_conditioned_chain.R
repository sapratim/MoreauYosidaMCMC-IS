
rm(list = ls())
source("nuclear_norm_functions.R")
iter_pc_chain <- 1e5
lamb_coeff <- 0.00001
sigma2_hat <- 0.01
alpha_hat <- 1.15/sigma2_hat
step_pcm <- 0.00008 

markov_chain <- mymala_cov_fn(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat,
                              iter = iter_pc_chain, delta = step_pcm)
start_value <- markov_chain[iter_pc_chain,]
save(start_value, file = "start_value.Rdata")
