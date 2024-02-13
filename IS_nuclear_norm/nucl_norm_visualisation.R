
rm(list = ls())
source("nuclear_norm_functions.R")
load("IS_MALA.Rdata")
load("Px_MALA.Rdata")
load("IS_HMC.Rdata")
load("Px_HMC.Rdata")
iter <- 1e5
##  IS samples
ismala_chain <- result_is[[1]]
ismala_weights <- result_is[[2]]
ishmc_chain <- result_ishmc[[1]]
ishmc_weights <- result_ishmc[[2]]

##  Px samples
pxmala_chain <- result_pxm
pxhmc_chain <- result_ishmc[[1]]

############# Posterior mean 

# post_mean <- function(samples, weights){
#   weight_mat <- matrix(0, nrow = iter, ncol = length(y))
#   for (i in 1:iter) {
#     weight_mat[i,] <- samples[i,]*exp(weights[i])
#   }
# num_sum <- apply(weight_mat, 2, sum)
# weights_sum <- sum(exp(weights))
# mean_val <- num_sum/weights_sum    #### Mean value
# return(mean_val)
# }
# 
# mean_ismala <- post_mean(ismala_chain, ismala_weights)
# mean_pxmala <- apply(pxmala_chain, 2 , mean)

################################ ESS evaluation MALA

mala_chain <- ismala_chain
weights <- ismala_weights
is_samp <- matrix(unlist(mala_chain), nrow = iter, ncol = length(y))
is_wts <- as.numeric(unlist(weights))
wts_mean <- mean(exp(is_wts))
asymp_var_is <- numeric(length = ncol(is_samp))
asymp_var_pxm <- numeric(length = ncol(is_samp))
for (j in 1:ncol(is_samp)) 
{
   num <- is_samp[, j]*exp(is_wts)
   sum_mat <- sum(num)
   is_est <- sum_mat / sum(exp(is_wts))
   input_mat <- cbind(num, exp(is_wts))  # input samples for mcse
   Sigma_mat <- mcse.multi(input_mat)$cov  # estimated covariance matrix of the tuple
   kappa_eta_mat <- cbind(1/wts_mean, is_est/wts_mean) # derivative of kappa matrix
   asymp_var_is[j] <- (kappa_eta_mat %*% Sigma_mat) %*% t(kappa_eta_mat) # IS asymptotic variance
}
#save(asymp_var_is, file = "asymp_var_is.Rdata")
comp_var <- mcse.mat(pxmala_chain)
asymp_var_pxm <- iter*(comp_var[,2]^2)   # PxMALA asymptotic variance
#save(asymp_var_pxm, file = "asymp_var_pxm.Rdata")
#load("asymp_var_is.Rdata")
rel_ess_mala <- asymp_var_pxm/asymp_var_is

################################ ESS evaluation HMC

hmc_chain <- ishmc_chain
weights <- ishmc_weights
is_samp <- matrix(unlist(hmc_chain), nrow = iter, ncol = length(y))
is_wts <- as.numeric(unlist(weights))
wts_mean <- mean(exp(is_wts))
asymp_var_ishmc <- numeric(length = ncol(is_samp))
asymp_var_pxhmc <- numeric(length = ncol(is_samp))
for (j in 1:ncol(is_samp)) 
{
  num <- is_samp[, j]*exp(is_wts)
  sum_mat <- sum(num)
  is_est <- sum_mat / sum(exp(is_wts))
  input_mat <- cbind(num, exp(is_wts))  # input samples for mcse
  Sigma_mat <- mcse.multi(input_mat)$cov  # estimated covariance matrix of the tuple
  kappa_eta_mat <- cbind(1/wts_mean, is_est/wts_mean) # derivative of kappa matrix
  asymp_var_ishmc[j] <- (kappa_eta_mat %*% Sigma_mat) %*% t(kappa_eta_mat) # IS asymptotic variance
}
#save(asymp_var_ishmc, file = "asymp_var_ishmc.Rdata")
comp_var <- mcse.mat(pxhmc_chain)
asymp_var_pxhmc <- iter*(comp_var[,2]^2)   # PxMALA asymptotic variance
#save(asymp_var_pxhmc, file = "asymp_var_pxhmc.Rdata")
rel_ess_hmc <- asymp_var_pxhmc/asymp_var_ishmc

pdf(file = "Hist_ess.pdf", width = 7, height = 4)
rel_ess_mala <- asymp_var_pxm/asymp_var_is
rel_ess_hmc <- asymp_var_pxhmc/asymp_var_ishmc
par(mfrow = c(1,2))
hist(rel_ess_mala)
hist(rel_ess_hmc)
dev.off()
