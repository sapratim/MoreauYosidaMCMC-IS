
source("nuclear_norm_functions.R")
load("start_value.Rdata")
iter <- 1e3
lamb_coeff <- 1e-5
sigma2_hat <- 0.01
alpha_hat <- 1.15/sigma2_hat
step_ismala <- 0.0001
step_pxmala <- 0.0001
step_isb <- 0.0001
step_pxb <- 0.0001
eps_is <- 0.0045
eps_px <-  0.005
L <- 10

parallel::detectCores()
num_cores <- 20
doParallel::registerDoParallel(cores = num_cores)
reps <- 100

output <- foreach(b = 1:reps) %dopar% {
######################### Covariance matrix MALA #########################
  
result_pxm <- px.mala(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat,
                        iter = iter, delta = step_pxmala, start = start_value)
comp_var <- mcse.mat(result_pxm)
asymp_var_pxm <- iter*(comp_var[,2]^2)   # PxMALA asymptotic variance
rm(result_pxm)

result_ism <- mymala(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat,
                      iter = iter, delta = step_ismala, start = start_value)
wts_mala <- exp(result_ism[[2]])
asymp_var_ism <- asymp_cov_func(result_ism[[1]], wts_mala)
n_eff_mala <- ((mean(wts_mala))^2/mean(wts_mala^2))
rm(result_ism)

######################### Covariance matrix Barker #########################

result_pxb <- px.barker(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat,
                      iter = iter, delta = step_pxb, start = start_value)
comp_var <- mcse.mat(result_pxb)
asymp_var_pxb <- iter*(comp_var[,2]^2)   # PxBarker asymptotic variance
rm(result_pxb)

result_isb <- mybarker(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat,
                    iter = iter, delta = step_isb, start = start_value)
wts_bark <- exp(result_isb[[2]])
asymp_var_isb <- asymp_cov_func(result_isb[[1]], wts_bark)
n_eff_bark <- ((mean(wts_bark))^2/mean(wts_bark^2))
rm(result_isb)

######################### Covariance matrix HMC #########################

result_pxhmc <- pxhmc(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat, 
                      iter = iter, eps_hmc = eps_px, L=L, start = start_value)
comp_var <- mcse.mat(result_pxhmc[[1]])
asymp_var_pxhmc <- iter*(comp_var[,2]^2)   # PxMALA asymptotic variance
rm(result_pxhmc)

result_ishmc <- myhmc(y=y, alpha = alpha_hat,lambda = lamb_coeff, sigma2 = sigma2_hat, 
                      iter = iter, eps_hmc = eps_is, L=L, start = start_value)
wts_hmc <- exp(result_ishmc[[2]])
asymp_var_ishmc <- asymp_cov_func(result_ishmc[[1]], wts_hmc)
n_eff_hmc <- ((mean(wts_hmc))^2/mean(wts_hmc^2))
rm(result_ishmc)

list(asymp_var_ism, asymp_var_pxm, asymp_var_isb, asymp_var_pxb,
                asymp_var_ishmc, asymp_var_pxhmc, n_eff_mala, n_eff_bark, n_eff_hmc)
}
save(output, file = "output_nucl_norm.Rdata")
