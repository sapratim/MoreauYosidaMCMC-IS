
source("IS_trendfiltering_functions.R")

iter <- 1e3
delta <- 0.00004
lamb_coeff <- 0.000001
D_mat <- getD(k=1, n=1e2, x)   #  D matrix
cormat <- mymala_cov_fn(y, alpha_hat, sigma2_hat, k=1, grid = x, iter = 1e4, delta = delta)[[2]]
delta_samp <- 0.000000032
mala.is <- mymala(y, alpha_hat, sigma2_hat, k=1, grid=x, iter, delta = delta_samp, cormat)
px_mala <- px.mala(y, alpha_hat, sigma2_hat, k=1, grid=x, iter, delta = delta_samp, cormat)

# Asymptotic variance

asymp_covmat_is <- matrix(0, length(y), length(y))
asymp_covmat_pxm <- matrix(0, length(y), length(y))

is_samp <- matrix(unlist(mala.is[1]), nrow = iter, ncol = length(y))
is_wts <- as.numeric(unlist(mala.is[2]))
wts_mean <- mean(is_wts)
num <- is_samp*is_wts
sum_mat <- apply(num, 2, sum)
is_est <- sum_mat / sum(is_wts)
input_mat <- cbind(num, is_wts)  # input samples for mcse
Sigma_mat <- mcse.multi(input_mat)$cov  # estimated covariance matrix of the tuple
kappa_eta_mat <- cbind(diag(1/wts_mean, length(y)), is_est/wts_mean) # derivative of kappa matrix
asymp_covmat_is <- (kappa_eta_mat %*% Sigma_mat) %*% t(kappa_eta_mat)

# PxMALA asymptotic variance

asymp_covmat_pxm <- mcse.multi(px_mala)$cov

# PxBarker asymptotic variance

# asymp_covmat_pxb[i] <- mcse.multi(px_bark)$cov

# Asymptotic variance comparison
var_mat <- cbind(det(asymp_covmat_is), det(asymp_covmat_pxm)) #asymp_covmat_pxb)
colnames(var_mat) <- c("Imp_sampling", "PxMala")

var_mat

