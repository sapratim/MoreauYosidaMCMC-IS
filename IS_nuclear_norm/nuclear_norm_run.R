rm(list = ls())
set.seed(8024248)
source("nuclear_norm_functions.R")
load("start_value.Rdata")
iter <- 1e5
lamb_coeff <- 1e-4
sigma2_hat <- 0.01
alpha_hat <- 1.15/sigma2_hat
step_ismala <- 0.00012
step_pxmala <- 0.0001
eps_is <- 0.008
eps_px <-  0.004
L <- 10

parallel::detectCores()
num_cores <- 10
doParallel::registerDoParallel(cores = num_cores)
reps <- 10

output <- foreach(b = 1:reps) %dopar% {
######################### ESS evaluation MALA #########################
  
result_pxm <- px.mala(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat,
                        iter = iter, delta = step_pxmala, start = start_value)
comp_var <- mcse.mat(result_pxm)
asymp_var_pxm <- iter*(comp_var[,2]^2)   # PxMALA asymptotic variance
rm(result_pxm)

result_is <- mymala(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat,
                      iter = iter, delta = step_ismala, start = start_value)
wts_mean <- mean(exp(result_is[[2]]))
asymp_var_is <- numeric(length = ncol(result_is[[1]]))
for (j in 1:ncol(result_is[[1]])) 
{
  num <- result_is[[1]][, j]*exp(result_is[[2]])
  sum_mat <- sum(num)
  is_est <- sum_mat / sum(exp(result_is[[2]]))
  input_mat <- cbind(num, exp(result_is[[2]]))  # input samples for mcse
  Sigma_mat <- mcse.multi(input_mat)$cov  # estimated covariance matrix of the tuple
  kappa_eta_mat <- cbind(1/wts_mean, -is_est/wts_mean) # derivative of kappa matrix
  asymp_var_is[j] <- (kappa_eta_mat %*% Sigma_mat) %*% t(kappa_eta_mat) # IS asymptotic variance
}
wts_mala <- exp(result_is[[2]])
n_eff_mala <- ((mean(wts_mala))^2/mean(wts_mala^2))
rm(result_is)

######################### ESS evaluation HMC #########################

result_pxhmc <- pxhmc(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat, 
                      iter = iter, eps_hmc = eps_px, L=L, start = start_value)
comp_var <- mcse.mat(result_pxhmc[[1]])
asymp_var_pxhmc <- iter*(comp_var[,2]^2)   # PxMALA asymptotic variance
rm(result_pxhmc)

result_ishmc <- myhmc(y=y, alpha = alpha_hat,lambda = lamb_coeff, sigma2 = sigma2_hat, 
                      iter = iter, eps_hmc = eps_is, L=L, start = start_value)
wts_mean <- mean(exp(result_ishmc[[2]]))
asymp_var_ishmc <- numeric(length = ncol(result_ishmc[[1]]))
for (j in 1:ncol(result_ishmc[[1]])) 
{
  num <- result_ishmc[[1]][, j]*exp(result_ishmc[[2]])
  sum_mat <- sum(num)
  is_est <- sum_mat / sum(exp(result_ishmc[[2]]))
  input_mat <- cbind(num, exp(result_ishmc[[2]]))  # input samples for mcse
  Sigma_mat <- mcse.multi(input_mat)$cov  # estimated covariance matrix of the tuple
  kappa_eta_mat <- cbind(1/wts_mean, -is_est/wts_mean) # derivative of kappa matrix
  asymp_var_ishmc[j] <- (kappa_eta_mat %*% Sigma_mat) %*% t(kappa_eta_mat) # IS asymptotic variance
}
wts_hmc <- exp(result_ishmc[[2]])
n_eff_hmc <- ((mean(wts_hmc))^2/mean(wts_hmc^2))
rm(result_ishmc)

list(asymp_var_is, asymp_var_pxm, asymp_var_ishmc, asymp_var_pxhmc, n_eff_mala, n_eff_hmc)
}
save(output, file = "output_nucl_norm.Rdata")
