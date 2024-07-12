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

##  Posterior mean
post_mean <- post_mean_fn(my.hmc[[1]], my.hmc[[2]])

#  Quantile visualisation
level <- 0.025
upper_quant <- quantile_func(hmc_chain, weights, level)[[1]]
lower_quant <- quantile_func(hmc_chain, weights, level)[[2]]
post_med <- quantile_func(hmc_chain, weights, level)[[3]]
list(post_mean, post_med, asymp_covmat_is, 
     asymp_covmat_pxhmc, upper_quant, lower_quant)
}
save(output_hmc, file = "output_hmc.Rdata")
