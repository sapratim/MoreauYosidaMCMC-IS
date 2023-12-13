
source("IS_trendf_functions_Pereyra.R")
load("covmat.Rdata")
load("MC_pcm.Rdata")

iter <- 1e4
lamb_coeff <- 0.0005
delta_samp <- .025
D_mat <- getD(k=1, n=1e2, x)   #  D matrix

mymala <- function(y, alpha, sigma2, k, grid, iter, delta, covmat)
{
  samp.mym <- matrix(0, nrow = iter, ncol = length(y))
  lambda <- lamb_coeff
  wts_is_est <- numeric(length = iter)
  beta_current <- markov_chain[1e4,]
    samp.mym[1,] <- beta_current
  g_val <- alpha*sum(abs(D_mat%*%beta_current)) + sum((y - beta_current)^2)/(2*sigma2)
  prox_val.curr <- prox_func(beta_current, lambda, alpha, sigma2, k, grid)
  g_lambda_val <- prox_arg(prox_val.curr, beta_current, lambda=lambda, y, sigma2, alpha)
  wts_is_est[1] <- g_lambda_val - g_val
  U <- sqrtm(covmat)
  mat.inv <- solve(covmat)
  accept <- 0
  for (i in 2:iter) 
  {
    beta_next <- beta_current +  ((delta / 2)*covmat)%*%log_gradpi(beta_current,lambda,y,sigma2,alpha,k,grid) + 
      (sqrt(delta)*U) %*% rnorm(length(beta_current), 0, 1)   # proposal step
    prox_val.next <- prox_func(beta_next, lambda, alpha, sigma2, k, grid)
    targ_val.next <- log_target(prox_val.next,beta_next,lambda,y,sigma2,alpha)
    targ_val.curr <- log_target(prox_val.curr,beta_current,lambda,y,sigma2,alpha)
    q.next_to_curr <- dmvnorm_fn(beta_current, beta_next + 
                                   ((delta / 2)*covmat)%*%log_gradpi(beta_next,lambda,y,sigma2,alpha,k,grid), mat.inv, delta)
    q.curr_to_next <- dmvnorm_fn(beta_next, beta_current + 
                                   ((delta / 2)*covmat)%*%log_gradpi(beta_current,lambda,y,sigma2,alpha,k,grid), mat.inv, delta) 
    mh.ratio <- targ_val.next + q.next_to_curr - (targ_val.curr + q.curr_to_next)  # mh  ratio
    # print(mh.ratio)
    if(log(runif(1)) <= mh.ratio)
    {
      samp.mym[i,] <- beta_next
      prox_val.curr <- prox_val.next
      g_val <- alpha*sum(abs(D_mat%*%beta_next)) + sum((y - beta_next)^2)/(2*sigma2)
      g_lambda_val <- prox_arg(prox_val.next, beta_next, lambda=lambda, y, sigma2, alpha)
      wts_is_est[i] <- g_lambda_val - g_val
      accept <- accept + 1
    }
    else
    {
      samp.mym[i,] <- beta_current
      g_val <- alpha*sum(abs(D_mat%*%beta_current)) + sum((y - beta_current)^2)/(2*sigma2)
      g_lambda_val <- prox_arg(prox_val.curr, beta_current, lambda=lambda, y, sigma2, alpha)
      wts_is_est[i] <- g_lambda_val - g_val
    }
    beta_current <- samp.mym[i,]
    if(i %% 10 == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(accept/iter)
  object <- list(samp.mym, wts_is_est)
  return(object)
}


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

