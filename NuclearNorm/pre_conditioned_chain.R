
rm(list = ls())
source("nuclear_norm_functions.R")
iter_pc_chain <- 1e5
lamb_coeff <- 0.00001
sigma2_hat <- 0.01
alpha_hat <- 1.15/sigma2_hat
step_pcm <- 0.00008 

markov_chain <- mymala_cov_fn(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat,
                              iter = iter_pc_chain, delta = step_pcm)

fn_val <- numeric(length = iter_pc_chain)
for (i in 1:iter_pc_chain) 
  {
  prox_val <- prox_func(markov_chain[i,], lamb_coeff, y, sigma2_hat, alpha_hat)
  fn_val[i] <- log_target(prox_val, markov_chain[i,], lamb_coeff,
                                         y, sigma2_hat, alpha_hat) 
}
index <- which.max(fn_val)
start_value <- markov_chain[index,]
save(start_value, file = "start_value.Rdata")
