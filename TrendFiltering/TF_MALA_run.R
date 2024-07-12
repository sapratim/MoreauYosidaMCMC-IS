
source("TF_functions.R")
load("pcm_last_iter.Rdata")

iter_mala <- 1e5
lamb_coeff <- 0.001
D_mat <- getD(k=1, n=1e2, x)   #  D matrix
delta_samp_is <- 0.0015
delta_samp_pxm <- 0.0008

output_mala <- list()
parallel::detectCores()
num_cores <- 50
doParallel::registerDoParallel(cores = num_cores)
reps <- 100

output <- foreach(b = 1:reps) %dopar% {
  mala.is <- mymala(y, alpha_hat, sigma2_hat, k=1, grid=x, iter = iter_mala, 
                    delta = delta_samp_is, start = pcm_last_iter)
  
  pxmala.run <- px.mala(y, alpha_hat, sigma2_hat, k=1, grid=x, iter = iter_mala, 
                        delta = delta_samp_pxm, start = pcm_last_iter)
 
mala_chain <- matrix(unlist(mala.is[[1]]), nrow = iter_mala, ncol = length(y))
weights <- exp(as.numeric(unlist(mala.is[[2]])))

# Asymptotic covariance matrix
asymp_covmat_is <- asymp_covmat_fn(mala_chain, weights) 
asymp_covmat_pxm <- mcse.multi(pxmala.run)$cov   

##  Posterior mean
post_mean <- post_mean_fn(mala.is[[1]], mala.is[[2]])

#  Quantile visualisation
level <- 0.025
upper_quant <- quantile_func(mala_chain, weights, level)[[1]]
lower_quant <- quantile_func(mala_chain, weights, level)[[2]]
post_med <- quantile_func(mala_chain, weights, level)[[3]]
list(post_mean, post_med, asymp_covmat_is, 
                asymp_covmat_pxm, upper_quant, lower_quant)
}
save(output_mala, file = "output_mala.Rdata")
