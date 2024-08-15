#####################################################################
################## Anisotropic Laplace output #######################
#####################################################################

source("ani_laplace_functions.R")
library(mcmcse)
library(foreach)
library(doParallel)

iter <- 1e6
reps <- 1e2
lamb_d5 <- seq(1e-2, 1, length = reps)
lamb_d10 <- seq(1e-2, 0.5, length = reps)
lamb_d20 <- seq(5e-3, 0.1, length = reps)
lamb_d50 <- seq(1e-3, 0.05, length = reps)

delta_is_d5 <- seq(0.11, 1.8, length = length(lamb_d5))
delta_px_d5 <- seq(0.11, 0.045, length = length(lamb_d5))

delta_is_d10 <- seq(0.027, 0.68, length = length(lamb_d10))
delta_px_d10 <- seq(0.022, 0.005, length = length(lamb_d10))  

delta_is_d20 <- seq(0.0075, 0.11, length = length(lamb_d20))
delta_px_d20 <- seq(0.0035, 0.0007, length = length(lamb_d20))  

delta_is_d50 <- seq(0.00095, 0.038, length = length(lamb_d50))
delta_px_d50 <- seq(0.00018, 0.00004, length = length(lamb_d50))  

output_laplace_d5 <- list()
output_laplace_d10 <- list()
output_laplace_d20 <- list()
output_laplace_d50 <- list()

parallel::detectCores()
num_cores <- 50
doParallel::registerDoParallel(cores = num_cores)

######################  dimension = 5  ################################

output_laplace_d5 <- foreach(lam_ani = 1:length(lamb_d5)) %dopar% {
  
  output_mala_d5 <- dimen_func(d <- 5, lambda <- lamb_d5[lam_ani], iter <- iter,
                               delta_is_d5[lam_ani], delta_px_d5[lam_ani])
  
  asymp_cov_ism_d5 <- asymp_covmat_fn(output_mala_d5[[1]][[1]], exp(output_mala_d5[[1]][[2]]))
  asymp_cov_pilam_d5 <- mcse.multi(output_mala_d5[[1]][[1]])$cov
  asymp_cov_pxm_d5 <- mcse.multi(output_mala_d5[[2]])$cov
  
  asymp_margvar_is5 <- diag(asymp_cov_ism_d5)   ### marginal variance ismala
  asymp_margvar_pilam5 <- diag(asymp_cov_pilam_d5)   ### marginal variance of is for pilambda
  rel_eff_5 <- diag(asymp_cov_pxm_d5)/ asymp_margvar_is5  ### relative efficiency wrt pi
  
  wts <- exp(output_mala_d5[[1]][[2]])
  ess_ism_d5 <- (mean(wts)^2)/mean(wts^2)
  
  list(asymp_margvar_is5, asymp_margvar_pilam5, rel_eff_5, ess_ism_d5)
}
save(output_laplace_d5, file = "output_d5.Rdata")

######################  dimension = 10  ################################

output_laplace_d10 <- foreach(lam_ani = 1:length(lamb_d10)) %dopar% {
  
  output_mala_d10 <- dimen_func(d <- 10, lambda <- lamb_d10[lam_ani], iter <- iter,
                                delta_is_d10[lam_ani], delta_px_d10[lam_ani])
  
  asymp_cov_ism_d10 <- asymp_covmat_fn(output_mala_d10[[1]][[1]], exp(output_mala_d10[[1]][[2]]))
  asymp_cov_pilam_d10 <- mcse.multi(output_mala_d10[[1]][[1]])$cov
  asymp_cov_pxm_d10 <- mcse.multi(output_mala_d10[[2]])$cov
  
  asymp_margvar_is10 <- diag(asymp_cov_ism_d10)   ### marginal variance ismala
  asymp_margvar_pilam10 <- diag(asymp_cov_pilam_d10)   ### marginal variance of is for pilambda
  rel_eff_10 <- diag(asymp_cov_pxm_d10)/ asymp_margvar_is10  ### relative efficiency wrt pi
  
  wts <- exp(output_mala_d10[[1]][[2]])
  ess_ism_d10 <- (mean(wts)^2)/mean(wts^2)
  
  list(asymp_margvar_is10, asymp_margvar_pilam10, rel_eff_10, ess_ism_d10)
}
save(output_laplace_d10, file = "output_d10.Rdata")

######################  dimension = 20  ################################

output_laplace_d20 <- foreach(lam_ani = 1:length(lamb_d20)) %dopar% {
  
  output_mala_d20 <- dimen_func(d <- 20, lambda <- lamb_d20[lam_ani], iter <- iter,
                                delta_is_d20[lam_ani], delta_px_d20[lam_ani])
  
  asymp_cov_ism_d20 <- asymp_covmat_fn(output_mala_d20[[1]][[1]], exp(output_mala_d20[[1]][[2]]))
  asymp_cov_pilam_d20 <- mcse.multi(output_mala_d20[[1]][[1]])$cov
  asymp_cov_pxm_d20 <- mcse.multi(output_mala_d20[[2]])$cov
  
  asymp_margvar_is20 <- diag(asymp_cov_ism_d20)   ### marginal variance ismala
  asymp_margvar_pilam20 <- diag(asymp_cov_pilam_d20)   ### marginal variance of is for pilambda
  rel_eff_20 <- diag(asymp_cov_pxm_d20)/ asymp_margvar_is20  ### relative efficiency wrt pi
  
  wts <- exp(output_mala_d20[[1]][[2]])
  ess_ism_d20 <- (mean(wts)^2)/mean(wts^2)
  
  list(asymp_margvar_is20, asymp_margvar_pilam20, rel_eff_20, ess_ism_d20)
}
save(output_laplace_d20, file = "output_d20.Rdata")

######################  dimension = 50  ################################

output_laplace_d50 <- foreach(lam_ani = 1:length(lamb_d50)) %dopar% {
  
  output_mala_d50 <- dimen_func(d <- 50, lambda <- lamb_d50[lam_ani], iter <- iter,
                                delta_is_d50[lam_ani], delta_px_d50[lam_ani])
  
  asymp_cov_ism_d50 <- asymp_covmat_fn(output_mala_d50[[1]][[1]], exp(output_mala_d50[[1]][[2]]))
  asymp_cov_pilam_d50 <- mcse.multi(output_mala_d50[[1]][[1]])$cov
  asymp_cov_pxm_d50 <- mcse.multi(output_mala_d50[[2]])$cov
  
  asymp_margvar_is50 <- diag(asymp_cov_ism_d50)   ### marginal variance ismala
  asymp_margvar_pilam50 <- diag(asymp_cov_pilam_d50)   ### marginal variance of is for pilambda
  rel_eff_50 <- diag(asymp_cov_pxm_d50)/ asymp_margvar_is50  ### relative efficiency wrt pi
  
  wts <- exp(output_mala_d50[[1]][[2]])
  ess_ism_d50 <- (mean(wts)^2)/mean(wts^2)
  
  list(asymp_margvar_is50, asymp_margvar_pilam50, rel_eff_50, ess_ism_d50)
}
save(output_laplace_d50, file = "output_d50.Rdata")
