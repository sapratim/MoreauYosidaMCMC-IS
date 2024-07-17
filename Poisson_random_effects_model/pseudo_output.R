
load("output_pseudo.Rdata")
library(mcmcse)
# to evaluate asymptotic covariance matrix
asymp_covmat_fn <- function(chain, weights)
{
  wts_mean <- mean(weights)
  num <- chain*weights
  sum_mat <- apply(num, 2, sum)
  is_est <- sum_mat / sum(weights)
  input_mat <- cbind(num, weights)  # input samples for mcse
  Sigma_mat <- mcse.multi(input_mat)$cov  # estimated covariance matrix of the tuple
  kappa_eta_mat <- cbind(diag(1/wts_mean, ncol(chain)), -is_est/wts_mean) # derivative of kappa matrix
  asymp_covmat <- (kappa_eta_mat %*% Sigma_mat) %*% t(kappa_eta_mat) # IS asymptotic variance
  return(asymp_covmat)
}

asymp_cov_ism <- asymp_covmat_fn(output_chain[[1]][[1]], exp(output_chain[[1]][[2]]))
asymp_cov_isb <- asymp_covmat_fn(output_chain[[3]][[1]], exp(output_chain[[3]][[2]]))
asymp_covmat_pxm <- mcse.multi(output_chain[[2]])$cov 
asymp_covmat_pxb <- mcse.multi(output_chain[[4]][[1]])$cov


#####  For univariate rel ess

# ratio of matrices
rel_ess_ratio_mala <- asymp_covmat_pxm/asymp_cov_ism
rel_ess_ratio_bark <- asymp_covmat_pxb/asymp_cov_isb

# actual univariate ess
uni_ess_mala <- diag(rel_ess_ratio_mala)
uni_ess_bark <- diag(rel_ess_ratio_bark)


# #  For posterior checks
# ism_mean <- apply(ismala[[1]], 2, mean)
# pxm_mean <- apply(pxmala, 2, mean)
# isb_mean <- apply(isbark[[1]], 2, mean)
# pxb_mean <- apply(pxbark[[1]], 2, mean)
