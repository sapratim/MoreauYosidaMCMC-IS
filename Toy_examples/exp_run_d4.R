#####################################################################
################# Exponential output visualisation #####################
#####################################################################

library(mcmcse)
library(foreach)
library(doParallel)
source("Toy_exp_functions.R")
curr_val <- 1

##########  opt lambda visualisation  ######################

num_cores <- 50
doParallel::registerDoParallel(cores = num_cores)

iter <- 1e6
reps <- 50

lamb_d4 <- seq(1e-4, 25, length = reps)

delta_is_d4 <- seq(1.2, 92, length = length(lamb_d4))
delta_px_d4 <- seq(1.4, 1.1, length = length(lamb_d4))


output_exp_d4 <- list()

output_exp_d4 <- foreach(lam_ani = 1:length(lamb_d4)) %dopar%{
  
  mym.output <- mymala(curr_val,iter, lamb_d4[lam_ani], delta_is_d4[lam_ani])
  pxm.output <- pxmala(curr_val, iter,lamb_d4[lam_ani],delta_px_d4[lam_ani])

  asymp_cov_ism_d4 <- asymp_covmat_fn(mym.output[[1]], exp(mym.output[[2]]))
  asymp_ess_pilam_d4 <- ess(mym.output[[1]], r = 1)
  asymp_cov_pxm_d4 <- mcse.multi(pxm.output, r = 1)$cov
  
  asymp_margvar_is4 <- diag(asymp_cov_ism_d4)   ### marginal variance ismala
  # asymp_margvar_pilam1 <- diag(asymp_cov_pilam_d4)   ### marginal variance of is for pilambda
  rel_eff_4 <- diag(asymp_cov_pxm_d4)/ asymp_margvar_is4  ### relative efficiency wrt pi
  
  wts <- exp(mym.output[[2]])
  ess_ism_d4 <- (mean(wts)^2)/mean(wts^2)
  
  list(asymp_margvar_is4, asymp_ess_pilam_d4/iter, rel_eff_4, ess_ism_d4)
}
save(output_exp_d4, file = "output_d4.Rdata")


# mym.output <- mymala(curr_val,iter, 1e-4, 1.2)
# pxm.output <- pxmala(curr_val, iter,25,1.1)
# 
# wts <- exp(mym.output[[2]])
# n_eff <- mean(wts)^2/mean(wts^2)
# n_eff
# 
