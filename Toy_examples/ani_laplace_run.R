#####################################################################
################## Anisotropic Laplace output #######################
#####################################################################

source("ani_laplace_functions.R")

iter <- 1e6
######################  dimension = 5  ################################
lamb_d5 <- seq(1e-2, 1, length = 1e2)
delta_is <- seq(0.025, 0.28, length = length(lamb_d5))
delta_px <- seq(0.025, 0.015, length = length(lamb_d5))

asymp_margvar_ism_d5 <- numeric(length = length(lamb_d5))
ess_ism_d5 <- numeric(length = length(lamb_d5))
asymp_var_pilam_d5 <- numeric(length = length(lamb_d5))

for (i in 1:length(lamb_d5)) 
  {
  output_mala_d5 <- dimen_func(d <- 5, lambda <- lamb_d5[i], iter <- iter, delta_is[i],
                               delta_px[i])
  
  asymp_cov_ism_5 <- asymp_covmat_fn(output_mala_d5[[1]][[1]], exp(output_mala_d5[[1]][[2]]))
  asymp_cov_pilam_d5 <- mcse.multi(output_mala_d5[[1]][[1]])$cov
  asymp_cov_pxm_5 <- mcse.multi(output_mala_d5[[2]])$cov
  
  asymp_margvar_ism_d5[i] <- mean(diag(asymp_cov_ism_5))   ####  Average asymptotic marginal variance
  asymp_var_pilam_d5[i] <- mean(diag(asymp_cov_pilam_d5))
  wts <- exp(output_mala_d5[[1]][[2]])
  ess_ism_d5[i] <- (mean(wts)^2)/mean(wts^2)
}

output_d5 <- list(asymp_margvar_ism_d5, ess_ism_d5, asymp_var_pilam_d5)

######################  dimension = 10  ##############################
lamb_d10 <- seq(1e-2, 0.5, length = 1e2)
delta_is <- seq(0.0025, 0.018, length = 1e2)
delta_px <- seq(0.0025, 0.0019, length = 1e2)  

asymp_margvar_ism_d10 <- numeric(length = length(lamb_d10))
ess_ism_d10 <- numeric(length = length(lamb_d10))
asymp_var_pilam_d10 <- numeric(length = length(lamb_d10))

for (i in 1:length(lamb_d10)) 
{
  output_mala_d10 <- dimen_func(d <- 10, lambda <- lamb_d10[i], iter <- iter, delta_is[i],
                               delta_px[i])
  
  asymp_cov_ism_10 <- asymp_covmat_fn(output_mala_d10[[1]][[1]], exp(output_mala_d10[[1]][[2]]))
  asymp_cov_pilam_d10 <- mcse.multi(output_mala_d10[[1]][[1]])$cov
  asymp_cov_pxm_10 <- mcse.multi(output_mala_d10[[2]])$cov
  
  asymp_margvar_ism_d10[i] <- mean(diag(asymp_cov_ism_10))   ####  Average asymptotic marginal variance
  asymp_var_pilam_d10[i] <- mean(diag(asymp_cov_pilam_d10))
  wts <- exp(output_mala_d10[[1]][[2]])
  ess_ism_d10[i] <- (mean(wts)^2)/mean(wts^2)
}

output_d10 <- list(asymp_margvar_ism_d10, ess_ism_d10, asymp_var_pilam_d10)

######################  dimension = 20  ##############################
lamb_d20 <- seq(5e-3, 0.1, length = 1e2)
delta_is <- seq(0.0003, 0.0007, length = 1e2)
delta_px <- seq(0.0003, 0.00025, length = 1e2)  

asymp_margvar_ism_d20 <- numeric(length = length(lamb_d20))
ess_ism_d20 <- numeric(length = length(lamb_d20))
asymp_var_pilam_d20 <- numeric(length = length(lamb_d20))

for (i in 1:length(lamb_d20)) 
{
  output_mala_d20 <- dimen_func(d <- 20, lambda <- lamb_d20[i], iter <- iter, delta_is[i],
                                delta_px[i])
  
  asymp_cov_ism_20 <- asymp_covmat_fn(output_mala_d20[[1]][[1]], exp(output_mala_d20[[1]][[2]]))
  asymp_cov_pilam_d20 <- mcse.multi(output_mala_d20[[1]][[1]])$cov
  asymp_cov_pxm_20 <- mcse.multi(output_mala_d20[[2]])$cov
  
  asymp_margvar_ism_d20[i] <- mean(diag(asymp_cov_ism_20))   ####  Average asymptotic marginal variance
  asymp_var_pilam_d20[i] <- mean(diag(asymp_cov_pilam_d20))
  wts <- exp(output_mala_d20[[1]][[2]])
  ess_ism_d20[i] <- (mean(wts)^2)/mean(wts^2)
}

output_d20 <- list(asymp_margvar_ism_d20, ess_ism_d20, asymp_var_pilam_d20)

######################  dimension = 50  ##############################
lamb_d50 <- seq(1e-3, 0.05, length = 1e2)
delta_is <- seq(0.00002, 0.00005, length = 1e2)
delta_px <- seq(0.00002, 0.000018, length = 1e2)  

asymp_margvar_ism_d50 <- numeric(length = length(lamb_d50))
ess_ism_d50 <- numeric(length = length(lamb_d50))
asymp_var_pilam_d50 <- numeric(length = length(lamb_d50))

for (i in 1:length(lamb_d50)) 
{
  output_mala_d50 <- dimen_func(d <- 50, lambda <- lamb_d50[i], iter <- iter, delta_is[i],
                                delta_px[i])
  
  asymp_cov_ism_50 <- asymp_covmat_fn(output_mala_d50[[1]][[1]], exp(output_mala_d50[[1]][[2]]))
  asymp_cov_pilam_d50 <- mcse.multi(output_mala_d50[[1]][[1]])$cov
  asymp_cov_pxm_50 <- mcse.multi(output_mala_d50[[2]])$cov
  
  asymp_margvar_ism_d50[i] <- mean(diag(asymp_cov_ism_50))   ####  Average asymptotic marginal variance
  rel_eff_matrix_50 <- asymp_cov_pxm_50/asymp_cov_ism_50
  asymp_var_pilam_d50[i] <- mean(diag(asymp_cov_pilam_d50))
  wts <- exp(output_mala_d50[[1]][[2]])
  ess_ism_d50[i] <- (mean(wts)^2)/mean(wts^2)
}
output_d50 <- list(asymp_margvar_ism_d50, ess_ism_d50, asymp_var_pilam_d50)

#######################################################################

output <- list(output_d5, output_d10, output_d20, output_d50)
save(output, file = "output_ani_laplace.Rdata")
