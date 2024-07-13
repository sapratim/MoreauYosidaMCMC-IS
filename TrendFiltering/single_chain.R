
source("TF_functions.R")
load("warmup.Rdata")

iter <- 1e5
lamb_coeff <- 0.001
delta_samp_is <- 0.0015
delta_samp_pxm <- 0.0008
delta_bark_my <- 0.0035
delta_bark_px <- 0.0025
D_mat <- getD(k=1, n=1e2, x)   #  D matrix
alpha_hat <- 5   # obtained from the first dataset
sigma2_hat <- 9  # obtained from the first dataset
k <- 1

wts_func <- function(beta,lambda,alpha,sigma2,k,grid)  # function for log of weights
{
  g_val <- alpha*sum(abs(D_mat%*%beta)) + sum((y - beta)^2)/(2*sigma2)
  proxval <- prox_func(beta, lambda, alpha, sigma2, k, grid)
  g_lambda_val <- prox_arg(proxval, beta, lambda=lambda, y, sigma2, alpha)
  wt_val <- g_lambda_val - g_val
  return(wt_val)
}

mala.is <- mymala(y, alpha_hat, sigma2_hat, k=1, grid=x, iter = iter, 
                  delta = delta_samp_is, start = warmup_end_iter)

pxmala.run <- px.mala(y, alpha_hat, sigma2_hat, k=1, grid=x, iter = iter, 
                      delta = delta_samp_pxm, start = warmup_end_iter)

mybark <- mybarker(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter_bark,
                   delta = delta_bark_my, start = warmup_end_iter)

pxbark <- px.barker(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter_bark,
                   delta = delta_bark_px, start = warmup_end_iter)

my.hmc <- myhmc(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter,
                eps_hmc = 0.015, L = 100, start = warmup_end_iter)

px.hmc <- pxhmc(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter,
                eps_hmc = 0.0003, L = 100, start = warmup_end_iter) 
output_single_run <- list(mala.is[[1]], pxmala.run, mybark[[1]], 
                              pxbark[[1]], my.hmc[[1]], px.hmc[[1]])

#######   Evaluation of log weights
wts_is_mala <- numeric(length = iter)
wts_is_bark <- numeric(length = iter)
wts_is_hmc <- numeric(length = iter)
for (i in 1:iter) {
  wts_is_mala[i] <- wts_func(output_single_run[[1]][i,], lamb_coeff, 
                             alpha_hat, sigma2_hat, k = k, grid = x)
  wts_is_bark[i] <- wts_func(output_single_run[[3]][i,], lamb_coeff, 
                             alpha_hat, sigma2_hat, k = k, grid = x)
  wts_is_hmc[i] <- wts_func(output_single_run[[5]][i,], lamb_coeff, 
                            alpha_hat, sigma2_hat, k = k, grid = x)
}
log_wts <- list(wts_is_mala, wts_is_bark, wts_is_hmc)
save(output_single_run, file = "single_chain.Rdata")
save(log_wts, file = "single_chain_log_weights.Rdata")
