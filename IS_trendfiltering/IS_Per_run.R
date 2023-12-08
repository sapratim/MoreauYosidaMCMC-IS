
source("IS_trendf_Pereyra.R")
load("covmat.Rdata")

iter <- 1e4
lamb_coeff <- 0.0005
delta_samp <- .042
D_mat <- getD(k=1, n=1e2, x)   #  D matrix

mala.is <- mymala(y, alpha_hat, sigma2_hat, k=1, grid=x, iter = iter, 
                                    delta = delta_samp, covmat = covmat)

mala_chain <- mala.is[[1]]
weights <- mala.is[[2]]
save(mala_chain, file = "mala_chain.Rdata")
save(weights, file = "weights.Rdata")

# Asymptotic variance

# asymp_covmat_is <- matrix(0, length(y), length(y))
# asymp_covmat_pxm <- matrix(0, length(y), length(y))
# 
# is_samp <- matrix(unlist(mala.is[1]), nrow = iter, ncol = length(y))
# is_wts <- as.numeric(unlist(mala.is[2]))
# wts_mean <- mean(is_wts)
# num <- is_samp*is_wts
# sum_mat <- apply(num, 2, sum)
# is_est <- sum_mat / sum(is_wts)
# input_mat <- cbind(num, is_wts)  # input samples for mcse
# Sigma_mat <- mcse.multi(input_mat)$cov  # estimated covariance matrix of the tuple
# kappa_eta_mat <- cbind(diag(1/wts_mean, length(y)), is_est/wts_mean) # derivative of kappa matrix
# asymp_covmat_is <- (kappa_eta_mat %*% Sigma_mat) %*% t(kappa_eta_mat)
# 
# # PxMALA asymptotic variance
# 
# asymp_covmat_pxm <- mcse.multi(px_mala)$cov
# 
# # PxBarker asymptotic variance
# 
# # asymp_covmat_pxb[i] <- mcse.multi(px_bark)$cov
# 
# # Asymptotic variance comparison
# var_mat <- cbind(sum(log(eigen(asymp_covmat_is)$values)), 
#                  sum(log(eigen(asymp_covmat_pxm)$values))) #asymp_covmat_pxb)
# colnames(var_mat) <- c("Imp_sampling", "PxMala")
# 
# var_mat

