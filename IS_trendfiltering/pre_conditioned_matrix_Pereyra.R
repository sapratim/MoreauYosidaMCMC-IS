source("IS_trendf_functions_Pereyra.R")

iter <- 1e4
lamb_coeff <- 0.0005
D_mat <- getD(k=1, n=1e2, x)   #  D matrix
delta <- .001

# MYMALA sampling for covariance matrix estimation

mymala_cov_fn <- function(y, alpha, sigma2, k, grid, iter, delta) #pre-conditioned mala
{
  samp.mym <- matrix(0, nrow = iter, ncol = length(y))
  lambda <- lamb_coeff
  beta_current <- y
  prox_val.curr <- prox_func_pcm(beta_current, lambda, alpha,sigma2, k, grid)
  samp.mym[1,] <- beta_current
  accept <- 0
  for (i in 2:iter) 
  {
    beta_next <- rnorm(length(beta_current), beta_current + 
                         (delta / 2)*log_gradpi(beta_current,lambda,y,sigma2,alpha,k,grid), 
                       sqrt(delta))   # proposal step
    prox_val.next <- prox_func_pcm(beta_next, lambda, alpha, sigma2, k, grid)
    targ_val.next <- log_target(prox_val.next,beta_next,lambda,y,sigma2,alpha)
    targ_val.curr <- log_target(prox_val.curr,beta_current,lambda,y,sigma2,alpha)
    q.next_to_curr <- sum(dnorm(beta_current, beta_next + 
                                  (delta / 2)*log_gradpi(beta_next,lambda,y,sigma2,alpha,k,grid),
                                sqrt(delta), log = TRUE))
    q.curr_to_next <- sum(dnorm(beta_next, beta_current + 
                                  (delta / 2)*log_gradpi(beta_current,lambda,y,sigma2,alpha,k,grid),
                                sqrt(delta), log = TRUE))
    mh.ratio <- targ_val.next + q.next_to_curr - (targ_val.curr + q.curr_to_next)  # mh  ratio
    # print(mh.ratio)
    if(log(runif(1)) <= mh.ratio)
    {
      samp.mym[i,] <- beta_next
      prox_val.curr <- prox_val.next
      accept <- accept + 1
    }
    else
    {
      samp.mym[i,] <- beta_current
    }
    beta_current <- samp.mym[i,]
    if(i %% 10 == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(accept/iter)
  object <- samp.mym
  cov_mat <- cov(object)
  result <- list(object, cov_mat)
  return(result)
}

prec_mat_chain <- mymala_cov_fn(y, alpha_hat, sigma2_hat, k=1, 
                                grid = x, iter, delta = delta)
markov_chain <- prec_mat_chain[[1]]
covmat <- prec_mat_chain[[2]]
save(markov_chain, file = "MC_pcm.Rdata")
save(covmat, file = "covmat.Rdata")

# px.mala <- function(y, alpha, sigma2, k, grid, iter, delta, covmat)
# {
#   samp.pxm <- matrix(0, nrow = iter, ncol = length(y))
#   lambda <- lamb_coeff*sigma2
#   beta_current <- trendfilter(grid,y, k=k,lambda = sigma2*alpha,
#                               control = trendfilter.control.list(obj_tol = tol, max_iter = 1e3L))$beta
#   samp.pxm[1,] <- beta_current
#   accept <- 0
#   U <- sqrtm(covmat)
#   mat.inv <- solve(covmat)
#   for (i in 2:iter) 
#   {
#     beta_next <- beta_current +  ((delta / 2)*covmat)%*%log_gradpi(beta_current,lambda,y,sigma2,alpha,k,grid) + 
#       (sqrt(delta)*U) %*% rnorm(length(beta_current), 0, 1)   # proposal step
#     U_betanext <- - (sum((y - beta_next)^2)/(2*sigma2) + alpha*(sum(abs(D_mat%*%beta_next))))
#     U_betacurr <- - (sum((y - beta_current)^2)/(2*sigma2) + alpha*(sum(abs(D_mat%*%beta_current))))
#     q.next_to_curr <- dmvnorm_fn(beta_current, beta_next + 
#                                    ((delta / 2)*covmat)%*%log_gradpi(beta_next,lambda,y,sigma2,alpha,k,grid), mat.inv, delta)
#     q.curr_to_next <- dmvnorm_fn(beta_next, beta_current + 
#                                    ((delta / 2)*covmat)%*%log_gradpi(beta_current,lambda,y,sigma2,alpha,k,grid), mat.inv, delta) 
#     mh.ratio <- U_betanext + q.next_to_curr - (U_betacurr + q.curr_to_next)
#     if(log(runif(1)) <= mh.ratio)
#     {
#       samp.pxm[i,] <- beta_next
#       accept <- accept + 1
#     }
#     else
#     {
#       samp.pxm[i,] <- beta_current
#     }
#     beta_current <- samp.pxm[i,]
#     if(i %% 1000 == 0){
#       print(i)
#     }
#   }
#   print(accept/iter)
#   return(samp.pxm)
# }


# px.barker <- function(in_val, iter, lambda, delta)
# {
#   samp.bark <- numeric(length = iter)
#   samp.bark[1] <- in_val
#   accept <- 0
#   for (i in 2:iter)
#   {
#     propval <- bark.prop(in_val, delta, lambda)
#     mh.ratio <- target_val(propval) + log(bark.dens(propval, in_val, delta, lambda)) - target_val(in_val) -
#       log(bark.dens(in_val, propval, delta, lambda))
#     if(log(runif(1)) <= mh.ratio)
#     {
#       samp.bark[i] <- propval
#       accept <- accept + 1
#     }
#     else
#     {
#       samp.bark[i] <- in_val
#     }
#     in_val <- samp.bark[i]
#   }
#   return(samp.bark)
# }

