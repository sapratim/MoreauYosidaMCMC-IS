
source("TF_functions.R")
load("pcm_last_iter.Rdata")

iter_bark <- 1e5
lamb_coeff <- 0.001
D_mat <- getD(k=1, n=1e2, x)   #  D matrix
delta_bark_my <- 0.06
delta_bark_px <- 0.042

output_bark <- list()

parallel::detectCores()
num_cores <- 50
doParallel::registerDoParallel(cores = num_cores)


output_bark <- foreach(b = 1:num_cores) %dopar% {
  pxbark <- px.barker(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter_bark,delta = delta_bark_px,
                      start = pcm_last_iter)
  mybark <- mybarker(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter_bark,delta = delta_bark_my,
                     start = pcm_last_iter)
  bark_chain <- matrix(unlist(mybark[[1]]), nrow = iter_bark, ncol = length(y))
  weights <- exp(as.numeric(unlist(mybark[[2]])))
  
  # Asymptotic covariance matrix
  asymp_covmat_is <- asymp_covmat_fn(bark_chain, weights) 
  asymp_covmat_pxb <- mcse.multi(pxbark[[1]])$cov   
  
  # Relative ESS
  rel_ess <- (det(asymp_covmat_pxb)/det(asymp_covmat_is))^(1/length(y))
  
##  Posterior mean
post_mean <- post_mean_fn(mybark[[1]],mybark[[2]])

#  Quantile visualisation

augm_mat <- cbind(mybark[[1]],mybark[[2]])

upper_quant <- numeric(length = length(y))
lower_quant <- numeric(length = length(y))
post_med <- numeric(length = length(y))
signif_level <- 0.025

for (i in 1:length(y)) 
{
  initial_mat <- quant(i, augm_mat)
  mat_sum <- apply(initial_mat, 2, sum)
  wts_prop <- initial_mat[,2]/mat_sum[2]
  final_mat <- cbind(initial_mat[1,], cumsum(wts_prop))
  lower_index <- min(which(final_mat[,2] >= signif_level))
  upper_index <- min(which(final_mat[,2] >= (1 - signif_level)))
  med_index <- min(which(final_mat[,2] >= 0.5))
  upper_quant[i] <- initial_mat[upper_index,1]
  lower_quant[i] <- initial_mat[lower_index,1]
  post_med[i] <- initial_mat[med_index,1]
}

acc_rate_is <- mybark[[3]]
acc_rate_pxb <- pxbark[[2]]
list(post_mean, post_med, asymp_covmat_is, asymp_covmat_pxbark, 
     upper_quant, lower_quant, rel_ess, acc_rate_is, acc_rate_pxb)
}

save(output_bark, file = "output_bark.Rdata")

