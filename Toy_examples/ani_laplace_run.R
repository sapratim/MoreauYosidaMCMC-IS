#####################################################################
################## Anisotropic Laplace output #######################
#####################################################################

source("ani_laplace_functions.R")
library(mcmcse)
library(foreach)
library(doParallel)

iter <- 5e3
output_laplace <- list()

parallel::detectCores()
num_cores <- 4
doParallel::registerDoParallel(cores = num_cores)
reps <- 1e2

lamb_d5 <- seq(1e-2, 1, length = reps)
lamb_d10 <- seq(1e-2, 0.5, length = reps)
lamb_d20 <- seq(5e-3, 0.1, length = reps)
lamb_d50 <- seq(1e-3, 0.05, length = reps)

#  For saving useful information viz. marginal variances and ess

asymp_margvar_is5 <- matrix(0, nrow = reps, ncol = 5)
asymp_margvar_pilam5 <- matrix(0, nrow = reps, ncol = 5)
rel_eff_5 <- matrix(0, nrow = reps, ncol = 5)

asymp_margvar_is10 <- matrix(0, nrow = reps, ncol = 10)
asymp_margvar_pilam10 <- matrix(0, nrow = reps, ncol = 10)
rel_eff_10 <- matrix(0, nrow = reps, ncol = 10)

asymp_margvar_is20 <- matrix(0, nrow = reps, ncol = 20)
asymp_margvar_pilam20 <- matrix(0, nrow = reps, ncol = 20)
rel_eff_20 <- matrix(0, nrow = reps, ncol = 20)

asymp_margvar_is50 <- matrix(0, nrow = reps, ncol = 50)
asymp_margvar_pilam50 <- matrix(0, nrow = reps, ncol = 50)
rel_eff_50 <- matrix(0, nrow = reps, ncol = 50)


##  For saving importance ess

ess_ism_d5 <- numeric(length = reps)
ess_ism_d10 <- numeric(length = reps)
ess_ism_d20 <- numeric(length = reps)
ess_ism_d50 <- numeric(length = reps)

delta_is_d5 <- seq(0.11, 1.8, length = length(lamb_d5))
delta_px_d5 <- seq(0.11, 0.045, length = length(lamb_d5))

delta_is_d10 <- seq(0.027, 0.68, length = length(lamb_d10))
delta_px_d10 <- seq(0.022, 0.005, length = length(lamb_d10))  

delta_is_d20 <- seq(0.0075, 0.11, length = length(lamb_d20))
delta_px_d20 <- seq(0.0035, 0.0007, length = length(lamb_d20))  

delta_is_d50 <- seq(0.00095, 0.038, length = length(lamb_d50))
delta_px_d50 <- seq(0.00018, 0.00004, length = length(lamb_d50))  

i <- 0
output_laplace <- foreach(b = 1:reps) %dopar% {
  ######################  dimension = 5  ################################
  i <- i+1
  
  output_mala_d5 <- dimen_func(d <- 5, lambda <- lamb_d5[i], iter <- iter, delta_is_d5[i],
                               delta_px_d5[i])
  asymp_cov_ism_d5 <- asymp_covmat_fn(output_mala_d5[[1]][[1]], exp(output_mala_d5[[1]][[2]]))
  asymp_cov_pilam_d5 <- mcse.multi(output_mala_d5[[1]][[1]])$cov
  asymp_cov_pxm_d5 <- mcse.multi(output_mala_d5[[2]])$cov

  asymp_margvar_is5[i,] <- diag(asymp_cov_ism_d5)
  asymp_margvar_pilam5[i,] <- diag(asymp_cov_pilam_d5)
  rel_eff_5[i,] <- diag(asymp_cov_pxm_d5)/ asymp_margvar_is5[i,]

  wts <- exp(output_mala_d5[[1]][[2]])
  ess_ism_d5[i] <- (mean(wts)^2)/mean(wts^2)
  
  ######################  dimension = 10  ################################
  
  output_mala_d10 <- dimen_func(d <- 10, lambda <- lamb_d10[i], iter <- iter, delta_is_d10[i],
                                delta_px_d10[i])
  
  asymp_cov_ism_d10 <- asymp_covmat_fn(output_mala_d10[[1]][[1]], exp(output_mala_d10[[1]][[2]]))
  asymp_cov_pilam_d10 <- mcse.multi(output_mala_d10[[1]][[1]])$cov
  asymp_cov_pxm_d10 <- mcse.multi(output_mala_d10[[2]])$cov
  
  asymp_margvar_is10[i,] <- diag(asymp_cov_ism_d10)
  asymp_margvar_pilam10[i,] <- diag(asymp_cov_pilam_d10)
  rel_eff_10[i,] <- diag(asymp_cov_pxm_d10)/ asymp_margvar_is10[i,]
  wts <- exp(output_mala_d10[[1]][[2]])
  ess_ism_d10[i] <- (mean(wts)^2)/mean(wts^2)
  
  ######################  dimension = 20  ################################
  
  output_mala_d20 <- dimen_func(d <- 20, lambda <- lamb_d20[i], iter <- iter, delta_is_d20[i],
                                delta_px_d20[i])
  
  asymp_cov_ism_d20 <- asymp_covmat_fn(output_mala_d20[[1]][[1]], exp(output_mala_d20[[1]][[2]]))
  asymp_cov_pilam_d20 <- mcse.multi(output_mala_d20[[1]][[1]])$cov
  asymp_cov_pxm_d20 <- mcse.multi(output_mala_d20[[2]])$cov
  
  asymp_margvar_is20[i,] <- diag(asymp_cov_ism_d20)
  asymp_margvar_pilam20[i,] <- diag(asymp_cov_pilam_d20)
  rel_eff_20[i,] <- diag(asymp_cov_pxm_d20)/ asymp_margvar_is20[i,]
  wts <- exp(output_mala_d20[[1]][[2]])
  ess_ism_d20[i] <- (mean(wts)^2)/mean(wts^2)
  
  ######################  dimension = 50  ################################
  
  output_mala_d50 <- dimen_func(d <- 50, lambda <- lamb_d50[i], iter <- iter, delta_is_d50[i],
                                delta_px_d50[i])
  
  asymp_cov_ism_d50 <- asymp_covmat_fn(output_mala_d50[[1]][[1]], exp(output_mala_d50[[1]][[2]]))
  asymp_cov_pilam_d50 <- mcse.multi(output_mala_d50[[1]][[1]])$cov
  asymp_cov_pxm_d50 <- mcse.multi(output_mala_d50[[2]])$cov
  
  asymp_margvar_is50[i,] <- diag(asymp_cov_ism_d50)
  asymp_margvar_pilam50[i,] <- diag(asymp_cov_pilam_d50)
  rel_eff_50[i,] <- diag(asymp_cov_pxm_d50)/ asymp_margvar_is50[i,]
  wts <- exp(output_mala_d50[[1]][[2]])
  ess_ism_d50[i] <- (mean(wts)^2)/mean(wts^2)
  
  list(asymp_margvar_is5, asymp_margvar_pilam5, rel_eff_5, ess_ism_d5, asymp_margvar_is10, 
       asymp_margvar_pilam10, rel_eff_10, ess_ism_d10, asymp_margvar_is20, asymp_margvar_pilam20,
       rel_eff_20, ess_ism_d20, asymp_margvar_is50, asymp_margvar_pilam50, rel_eff_50, ess_ism_d50)
}

save(output_laplace, file = "output_laplace.Rdata")
