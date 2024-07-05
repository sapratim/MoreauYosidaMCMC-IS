source("IS_trendf_functions_Pereyra.R")
load("single_chain.Rdata")

chain_length <- 1e5
lamb_coeff <- 0.001
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

#######   Evaluation of log weights
wts_is_mala <- numeric(length = chain_length)
wts_is_hmc <- numeric(length = chain_length)
for (i in 1:chain_length) {
  wts_is_mala[i] <- wts_func(output_single_run[[1]][i,], lamb_coeff, 
                             alpha_hat, sigma2_hat, k = k, grid = x)
  wts_is_hmc[i] <- wts_func(output_single_run[[3]][i,], lamb_coeff, 
                            alpha_hat, sigma2_hat, k = k, grid = x)
}

log_wts <- list(wts_is_mala, wts_is_hmc)

save(log_wts, file = "single_chain_log_weights.Rdata")
