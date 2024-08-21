#####################################################################
################# Exponential output visualisation ##################
#####################################################################
set.seed(100)
source("exp_d4_functions.R")
library(mcmcse)
library(foreach)
library(doParallel)
num_cores <- 50
doParallel::registerDoParallel(cores = num_cores)

iter <- 1e6
reps <- 50

lamb_d1 <- seq(1e-4, 5, length = reps)
lamb_d5 <- seq(1e-4, .6, length = reps)
lamb_d10 <- seq(1e-3, 0.5, length = reps)
lamb_d20 <- seq(1e-3, 0.2, length = reps)
# lamb_d50 <- seq(1e-2, 1.2, length = reps)

delta_is_d1 <- seq(.7, 14, length = length(lamb_d1))
delta_px_d1 <- seq(.7, .45, length = length(lamb_d1))

delta_is_d5 <- seq(0.35, 1.9, length = length(lamb_d5))
delta_px_d5 <- seq(0.35, 0.17, length = length(lamb_d5))

delta_is_d10 <- seq(0.2, 1, length = length(lamb_d10))
delta_px_d10 <- seq(0.2, 0.07, length = length(lamb_d10))  

delta_is_d20 <- seq(0.14, 0.45, length = length(lamb_d20))
delta_px_d20 <- seq(0.14, 0.04, length = length(lamb_d20))  

# delta_is_d50 <- seq(0.23, 1, length = length(lamb_d50))
# delta_px_d50 <- seq(0.23, 0.16, length = length(lamb_d50))  

output_exp_d1 <- list()
output_exp_d5 <- list()
output_exp_d10 <- list()
output_exp_d20 <- list()
# output_exp_d50 <- list()


######################  dimension = 1  ################################

output_exp_d1 <- foreach(lam_ani = 1:length(lamb_d1)) %dopar%{
  
  output_mala_d1 <- dimen_func(d = 1, lambda = lamb_d1[lam_ani], iter = iter,
                               delta_is_d1[lam_ani], delta_px_d1[lam_ani])
  
  asymp_cov_ism_d1 <- asymp_covmat_fn(output_mala_d1[[1]][[1]], exp(output_mala_d1[[1]][[2]]))
  asymp_ess_pilam_d1 <- ess(output_mala_d1[[1]][[1]], r = 1)
  asymp_cov_pxm_d1 <- mcse.multi(output_mala_d1[[2]], r = 1)$cov
  
  asymp_margvar_is1 <- diag(asymp_cov_ism_d1)   ### marginal variance ismala
  # asymp_margvar_pilam1 <- diag(asymp_cov_pilam_d1)   ### marginal variance of is for pilambda
  rel_eff_1 <- diag(asymp_cov_pxm_d1)/ asymp_margvar_is1  ### relative efficiency wrt pi
  
  wts <- exp(output_mala_d1[[1]][[2]])
  ess_ism_d1 <- (mean(wts)^2)/mean(wts^2)
  
  list(asymp_margvar_is1, asymp_ess_pilam_d1/iter, rel_eff_1, ess_ism_d1)
}
save(output_exp_d1, file = "output_d1.Rdata")


######################  dimension = 5  ################################

output_exp_d5 <- foreach(lam_ani = 1:length(lamb_d5)) %dopar%{
  
  output_mala_d5 <- dimen_func(d = 5, lambda = lamb_d5[lam_ani], iter = iter,
                               delta_is_d5[lam_ani], delta_px_d5[lam_ani])
  
  asymp_cov_ism_d5 <- asymp_covmat_fn(output_mala_d5[[1]][[1]], exp(output_mala_d5[[1]][[2]]))
  asymp_ess_pilam_d5 <- ess(output_mala_d5[[1]][[1]], r = 1)
  asymp_cov_pxm_d5 <- mcse.multi(output_mala_d5[[2]], r = 1)$cov
  
  asymp_margvar_is5 <- diag(asymp_cov_ism_d5)   ### marginal variance ismala
  # asymp_margvar_pilam5 <- diag(asymp_cov_pilam_d5)   ### marginal variance of is for pilambda
  rel_eff_5 <- diag(asymp_cov_pxm_d5)/ asymp_margvar_is5  ### relative efficiency wrt pi
  
  wts <- exp(output_mala_d5[[1]][[2]])
  ess_ism_d5 <- (mean(wts)^2)/mean(wts^2)
  
  list(asymp_margvar_is5, asymp_ess_pilam_d5/iter , rel_eff_5, ess_ism_d5)
}
save(output_exp_d5, file = "output_d5.Rdata")

######################  dimension = 10  ################################

output_exp_d10 <- foreach(lam_ani = 1:length(lamb_d10)) %dopar% {
  
  output_mala_d10 <- dimen_func(d = 10, lambda = lamb_d10[lam_ani], iter = iter,
                                delta_is_d10[lam_ani], delta_px_d10[lam_ani])
  
  asymp_cov_ism_d10 <- asymp_covmat_fn(output_mala_d10[[1]][[1]], exp(output_mala_d10[[1]][[2]]))
  asymp_ess_pilam_d10 <- ess(output_mala_d10[[1]][[1]], r = 1)
  asymp_cov_pxm_d10 <- mcse.multi(output_mala_d10[[2]], r = 1)$cov
  
  asymp_margvar_is10 <- diag(asymp_cov_ism_d10)   ### marginal variance ismala
  #asymp_margvar_pilam10 <- diag(asymp_cov_pilam_d10)   ### marginal variance of is for pilambda
  rel_eff_10 <- diag(asymp_cov_pxm_d10)/ asymp_margvar_is10  ### relative efficiency wrt pi
  
  wts <- exp(output_mala_d10[[1]][[2]])
  ess_ism_d10 <- (mean(wts)^2)/mean(wts^2)
  
  list(asymp_margvar_is10, asymp_ess_pilam_d10/iter , rel_eff_10, ess_ism_d10)
}
save(output_exp_d10, file = "output_d10.Rdata")

######################  dimension = 20  ################################

output_exp_d20 <- foreach(lam_ani = 1:length(lamb_d20)) %dopar% {
  
  output_mala_d20 <- dimen_func(d = 20, lambda = lamb_d20[lam_ani], iter = iter,
                                delta_is_d20[lam_ani], delta_px_d20[lam_ani])
  
  asymp_cov_ism_d20 <- asymp_covmat_fn(output_mala_d20[[1]][[1]], exp(output_mala_d20[[1]][[2]]))
  asymp_ess_pilam_d20 <- ess(output_mala_d20[[1]][[1]], r = 1)
  asymp_cov_pxm_d20 <- mcse.multi(output_mala_d20[[2]], r = 1)$cov
  
  asymp_margvar_is20 <- diag(asymp_cov_ism_d20)   ### marginal variance ismala
  #  asymp_margvar_pilam20 <- diag(asymp_cov_pilam_d20)   ### marginal variance of is for pilambda
  rel_eff_20 <- diag(asymp_cov_pxm_d20)/ asymp_margvar_is20  ### relative efficiency wrt pi
  
  wts <- exp(output_mala_d20[[1]][[2]])
  ess_ism_d20 <- (mean(wts)^2)/mean(wts^2)
  
  list(asymp_margvar_is20, asymp_ess_pilam_d20/iter, rel_eff_20, ess_ism_d20)
}
save(output_exp_d20, file = "output_d20.Rdata")

######################  dimension = 50  ################################

# output_exp_d50 <- foreach(lam_ani = 1:length(lamb_d50)) %dopar% {

#   output_mala_d50 <- dimen_func(d = 50, lambda = lamb_d50[lam_ani], iter = iter,
#                                 delta_is_d50[lam_ani], delta_px_d50[lam_ani])

#   asymp_cov_ism_d50 <- asymp_covmat_fn(output_mala_d50[[1]][[1]], exp(output_mala_d50[[1]][[2]]))
#   asymp_ess_pilam_d50 <- ess(output_mala_d50[[1]][[1]], r = 1)
#   asymp_cov_pxm_d50 <- mcse.multi(output_mala_d50[[2]], r = 1)$cov

#   asymp_margvar_is50 <- diag(asymp_cov_ism_d50)   ### marginal variance ismala
#  # asymp_margvar_pilam50 <- diag(asymp_cov_pilam_d50)   ### marginal variance of is for pilambda
#   rel_eff_50 <- diag(asymp_cov_pxm_d50)/ asymp_margvar_is50  ### relative efficiency wrt pi

#   wts <- exp(output_mala_d50[[1]][[2]])
#   ess_ism_d50 <- (mean(wts)^2)/mean(wts^2)

#   list(asymp_margvar_is50, asymp_ess_pilam_d50, rel_eff_50, ess_ism_d50)
# }
# save(output_exp_d50, file = "output_d50.Rdata")
