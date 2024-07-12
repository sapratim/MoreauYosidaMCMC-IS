
source("TF_functions.R")
load("pcm_last_iter.Rdata")

iter_mala <- 1e5
lamb_coeff <- 0.001
D_mat <- getD(k=1, n=1e2, x)   #  D matrix
delta_samp_is <- 0.00015
delta_samp_pxm <- 0.00012

output <- list()

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

# Relative ESS
rel_ess <- (det(asymp_covmat_pxm)/det(asymp_covmat_is))^(1/length(y))

##  Posterior mean
post_mean <- post_mean_fn(mala.is[[1]], mala.is[[2]])

#  Quantile visualisation

augm_mat <- cbind(mala.is[[1]],mala.is[[2]])

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

list(post_mean, post_med, Sigma_mat, 
     asymp_covmat_is, asymp_covmat_pxm, upper_quant, lower_quant, rel_ess)
}

save(output, file = "output_mala.Rdata")

# save(mala_chain, file = "mala_chain.Rdata")
# save(weights, file = "weights.Rdata")
# save(pxmala_chain, file = "pxmala_chain.Rdata")
