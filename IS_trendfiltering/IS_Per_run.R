
source("IS_trendf_Pereyra.R")

iter <- 1e4

lamb_coeff <- 0.0005
D_mat <- getD(k=1, n=1e2, x)   #  D matrix
delta <- .0006
prec_mat_chain <- mymala_cov_fn(y, alpha_hat, sigma2_hat, k=1, grid = x, iter, delta = delta)
covmat <- prec_mat_chain[[2]]
delta_samp <- .07
mala.is <- mymala(y, alpha_hat, sigma2_hat, k=1, grid=x, iter, delta = delta_samp, covmat)

freq_mode <- trendfilter(x,y, k=1,lambda = sigma2_hat*alpha_hat,
                  control = trendfilter.control.list(obj_tol = tol, max_iter = 1e3L))$beta
proxval <- prox_func(freq_mode, lamb_coeff,alpha_hat, sigma2_hat, k=1, grid = x)
mode_diff <- abs(freq_mode-proxval)
mode_diff   #   difference in modes

i <- 1

# Repeated execution gives density functions for different components
plot(density(mala.is[[1]][, i]))
abline(v = proxval[i], col = "red")
abline(v= freq_mode[i], col = "blue")
legend("topright", c("est_density_MCMC", "frequentist mode", "prox value at mode"), lty = 1,
       col = c("black", "blue", "red"), cex = 0.6, bty = "n")
i <- i + 1


###  Trace plots
plot.ts(mala.is[[1]][, 1:10])
plot.ts(mala.is[[1]][, 11:20])
plot.ts(mala.is[[1]][, 21:30])
plot.ts(mala.is[[1]][, 31:40])
plot.ts(mala.is[[1]][, 41:50])
plot.ts(mala.is[[1]][, 51:60])
plot.ts(mala.is[[1]][, 61:70])
plot.ts(mala.is[[1]][, 71:80])
plot.ts(mala.is[[1]][, 81:90])
plot.ts(mala.is[[1]][, 91:100])

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

