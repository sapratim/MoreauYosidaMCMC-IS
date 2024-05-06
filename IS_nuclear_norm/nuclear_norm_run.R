rm(list = ls())
source("nuclear_norm_functions.R")
load("start_value.Rdata")
iter <- 1e5 
lamb_coeff <- 1e-6
sigma2_hat <- 0.01
alpha_hat <- 1.15/sigma2_hat
step_ismala <- 0.00008
step_pxmala <- 0.00008

parallel::detectCores()
num_cores <- 50
doParallel::registerDoParallel(cores = num_cores)
reps <- 50

output <- foreach(b = 1:reps) %dopar% {
######################### ESS evaluation MALA #########################
result_is <- mymala(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat,
                      iter = iter, delta = step_ismala, start = start_value)
result_pxm <- px.mala(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat,
                      iter = iter, delta = step_pxmala, start = start_value)
is_samp <- matrix(unlist(result_is[[1]]), nrow = iter, ncol = length(y))
is_wts <- as.numeric(unlist(result_is[[2]]))
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
  kappa_eta_mat <- cbind(1/wts_mean, -is_est/wts_mean) # derivative of kappa matrix
  asymp_var_is[j] <- (kappa_eta_mat %*% Sigma_mat) %*% t(kappa_eta_mat) # IS asymptotic variance
}

   comp_var <- mcse.mat(result_pxm)
   asymp_var_pxm <- iter*(comp_var[,2]^2)   # PxMALA asymptotic variance
   rel_ess_mala <- asymp_var_pxm/asymp_var_is

######################### ESS evaluation HMC #########################

result_ishmc <- myhmc(y=y, alpha = alpha_hat,lambda = lamb_coeff, sigma2 = sigma2_hat, 
                      iter = iter, eps_hmc = 0.004, L=10, start = start_value)
result_pxhmc <- pxhmc(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat, 
                      iter = iter, eps_hmc = 0.004, L=10, start = start_value)
is_samp <- matrix(unlist(result_ishmc[[1]]), nrow = iter, ncol = length(y))
is_wts <- as.numeric(unlist(result_ishmc[[2]]))
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
  kappa_eta_mat <- cbind(1/wts_mean, -is_est/wts_mean) # derivative of kappa matrix
  asymp_var_ishmc[j] <- (kappa_eta_mat %*% Sigma_mat) %*% t(kappa_eta_mat) # IS asymptotic variance
}

    comp_var <- mcse.mat(result_pxhmc[[1]])
    asymp_var_pxhmc <- iter*(comp_var[,2]^2)   # PxMALA asymptotic variance
    rel_ess_hmc <- asymp_var_pxhmc/asymp_var_ishmc

  list(rel_ess_mala, rel_ess_hmc)
}
save(output, file = "output_nucl_norm.Rdata")
