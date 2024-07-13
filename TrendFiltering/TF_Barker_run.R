
source("TF_functions.R")
load("warmup.Rdata")

iter_bark <- 1e5
lamb_coeff <- 0.001
D_mat <- getD(k=1, n=1e2, x)   #  D matrix
delta_bark_my <- 0.0035
delta_bark_px <- 0.0025

output_bark <- list()

parallel::detectCores()
num_cores <- 50
doParallel::registerDoParallel(cores = num_cores)


output_bark <- foreach(b = 1:num_cores) %dopar% {
  mybark <- mybarker(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter_bark,delta = delta_bark_my,
                     start = warmup_end_iter)
  pxbark <- px.barker(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter_bark,delta = delta_bark_px,
                      start = warmup_end_iter)
  bark_chain <- matrix(unlist(mybark[[1]]), nrow = iter_bark, ncol = length(y))
  weights <- exp(as.numeric(unlist(mybark[[2]])))
  imp_ess <- (mean(weights)^2)/mean(weights^2)
  
# Asymptotic covariance matrix
asymp_covmat_is <- asymp_covmat_fn(bark_chain, weights) 
asymp_covmat_pxb <- mcse.multi(pxbark[[1]])$cov   
 
#  Posterior mean
post_mean <- post_mean_fn(mybark[[1]],mybark[[2]])

#  Quantile visualisation
level <- 0.025
upper_quant <- quantile_func(bark_chain, weights, level)[[1]]
lower_quant <- quantile_func(bark_chain, weights, level)[[2]]
post_med <- quantile_func(bark_chain, weights, level)[[3]]
list(post_mean, post_med, asymp_covmat_is, 
     asymp_covmat_pxb, upper_quant, lower_quant, imp_ess)
}
save(output_bark, file = "output_bark.Rdata")
