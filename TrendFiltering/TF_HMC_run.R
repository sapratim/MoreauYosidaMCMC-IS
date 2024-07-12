source("TF_functions.R")
load("pcm_last_iter.Rdata")

iter_hmc <- 1e5
lamb_coeff <- 0.001
D_mat <- getD(k=1, n=1e2, x)   #  D matrix

output_hmc <- list()

parallel::detectCores()
num_cores <- 50
doParallel::registerDoParallel(cores = num_cores)
reps <- 100

output_hmc <- foreach(b = 1:reps) %dopar% {
  my.hmc <- myhmc(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter_hmc,
                  eps_hmc = 0.015, L = 100, start = pcm_last_iter)
  
  px.hmc <- pxhmc(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter_hmc,
                  eps_hmc = 0.0003, L = 100, start = pcm_last_iter) 
hmc_chain <- matrix(unlist(my.hmc[[1]]), nrow = iter_hmc, ncol = length(y))
weights <- as.numeric(unlist(my.hmc[[2]]))

# Asymptotic covariance matrix
asymp_covmat_is <- asymp_covmat_fn(hmc_chain, weights) 
asymp_covmat_pxhmc <- mcse.multi(px.hmc[[1]])$cov   # PxMALA asymptotic variance

# Relative ESS
rel_ess <- (det(asymp_covmat_pxhmc)/det(asymp_covmat_is))^(1/length(y))

##  Posterior mean
post_mean <- post_mean_fn(my.hmc[[1]], my.hmc[[2]])

#  Quantile visualisation

augm_mat <- cbind(my.hmc[[1]],my.hmc[[2]])

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

acc_rate_is <- my.hmc[[3]]
acc_rate_pxhmc <- px.hmc[[2]]
list(post_mean, post_med, Sigma_mat, asymp_covmat_is, asymp_covmat_pxhmc, 
     upper_quant, lower_quant, rel_ess, acc_rate_is, acc_rate_pxhmc)
}

save(output_hmc, file = "output_hmc.Rdata")
