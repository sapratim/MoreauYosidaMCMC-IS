#####################################################################
########### Anisotropic Laplace output visualisation ################
#####################################################################

source("ani_laplace_functions.R")
library(mcmcse)

lambda <- 1e-5
iter <- 1e4
######################  dimension = 5  ##############################
d <- 5
beta <- seq(1, d, by = 1)
start <- rep(0.5, d)

mym.output <- mymala(beta, start, lambda, iter, delta <- 0.003)
pxm.output <- pxmala(beta, start, lambda, iter, delta <- 0.0025)

asymp_cov_ism_5 <- asymp_covmat_fn(mym.output[[1]], exp(mym.output[[2]]))
asymp_cov_pxm_5 <- mcse.multi(pxm.output)$cov

rel_eff_matrix_5 <- asymp_cov_pxm_5/asymp_cov_ism_5
avg_rel_eff_5 <- mean(diag(rel_eff_matrix_5))

######################  dimension = 10  ##############################
d <- 10
beta <- seq(1, d, by = 1)
start <- rep(0.5, d)

mym.output <- mymala(beta, start, lambda, iter, delta <- 0.0025)
pxm.output <- pxmala(beta, start, lambda, iter, delta <- 0.0025)

asymp_cov_ism_10 <- asymp_covmat_fn(mym.output[[1]], exp(mym.output[[2]]))
asymp_cov_pxm_10 <- mcse.multi(pxm.output)$cov

rel_eff_matrix_10 <- asymp_cov_pxm_10/asymp_cov_ism_10
avg_rel_eff_10 <- mean(diag(rel_eff_matrix_10))

######################  dimension = 20  ##############################
d <- 20
beta <- seq(1, d, by = 1)
start <- rep(0.5, d)

mym.output <- mymala(beta, start, lambda, iter, delta <- 0.00025)
pxm.output <- pxmala(beta, start, lambda, iter, delta <- 0.00025)

asymp_cov_ism_20 <- asymp_covmat_fn(mym.output[[1]], exp(mym.output[[2]]))
asymp_cov_pxm_20 <- mcse.multi(pxm.output)$cov

rel_eff_matrix_20 <- asymp_cov_pxm_20/asymp_cov_ism_20
avg_rel_eff_20 <- mean(diag(rel_eff_matrix_20))

######################  dimension = 50  ##############################
d <- 50
beta <- seq(1, d, by = 1)
start <- rep(0.5, d)

mym.output <- mymala(beta, start, lambda, iter, delta <- 0.00005)
pxm.output <- pxmala(beta, start, lambda, iter, delta <- 0.00005)

asymp_cov_ism_50 <- asymp_covmat_fn(mym.output[[1]], exp(mym.output[[2]]))
asymp_cov_pxm_50 <- mcse.multi(pxm.output)$cov

rel_eff_matrix_50 <- asymp_cov_pxm_50/asymp_cov_ism_50
avg_rel_eff_50 <- mean(diag(rel_eff_matrix_50))
