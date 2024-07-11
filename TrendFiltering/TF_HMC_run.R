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

hmc_chain <- my.hmc[[1]]
weights <- my.hmc[[2]]
pxhmc_chain <- px.hmc[[1]]
is_samp <- matrix(unlist(hmc_chain), nrow = iter_hmc, ncol = length(y))
is_wts <- as.numeric(unlist(weights))
wts_mean <- mean(exp(is_wts))
num <- is_samp*exp(is_wts)
sum_mat <- apply(num, 2, sum)
is_est <- sum_mat / sum(exp(is_wts))
input_mat <- cbind(num, exp(is_wts))  # input samples for mcse
Sigma_mat <- mcse.multi(input_mat)$cov  # estimated covariance matrix of the tuple
kappa_eta_mat <- cbind(diag(1/wts_mean, length(y)), -is_est/wts_mean) # derivative of kappa matrix

asymp_covmat_is <- (kappa_eta_mat %*% Sigma_mat) %*% t(kappa_eta_mat) # IS asymptotic variance

asymp_covmat_pxhmc <- mcse.multi(pxhmc_chain)$cov   # PxMALA asymptotic variance

rel_ess <- (det(asymp_covmat_pxhmc)/det(asymp_covmat_is))^(1/length(y))

##  Posterior mean

weight_mat <- matrix(0, nrow = iter_hmc, ncol = length(y))
for (i in 1:iter_hmc) {
  weight_mat[i,] <- hmc_chain[i,]*exp(weights[i])
}
num_sum <- apply(weight_mat, 2, sum)
weights_sum <- sum(exp(weights))
post_mean <- num_sum/weights_sum

#  Quantile visualisation

augm_mat <- cbind(hmc_chain,weights)

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
