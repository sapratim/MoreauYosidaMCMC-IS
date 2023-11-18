
## Code for asymptotic variance comparison of weighted importance sampling estimator 
## for trend filtering example
library(mcmcse)
library(coda)
library(glmgen)
library(Matrix)
library(spmrf)
library(expm)

set.seed(301297)
alpha_hat <- 5   # obtained from the first dataset
sigma2_hat <- 12.5  # obtained from the first dataset
x <- seq(1,100,len=100)
f <- Vectorize(function(x){if(x<=35){x} else if(x<=70){70-x} else{0.5*x-35}})
fx_linear <- f(x)
y <- fx_linear + rnorm(length(x), sd = 3)

# function calculates the inside of the proximal function

getD <- function(k, n, x=NULL){
  if(is.null(x)){
    x <- 1:n
  }
  diags <- list(rep(-1,n),rep(1,n))
  D <- Matrix::bandSparse(n-1,n,k=c(0,1),diag=diags,symm=F)
  if(k>=1){
    for(i in 1:k){
      leftD <- Matrix::bandSparse(n-i-1,n-i,k=c(0,1),diag=diags,symm=F)
      xdiag <- Matrix::Diagonal(n-i,i/diff(x,lag=i))
      D <- leftD %*% xdiag %*% D
    }
  }
  return(D)
}

log_target <- function(eta,beta,lambda,y,sigma2,alpha)
{
  dens_val <- alpha*(sum(abs(D_mat%*%eta))) + sum((beta-eta)^2)/(2*lambda) + 
                             sum((y-eta)^2)/(2*sigma2)
  return(-dens_val)
}

prox_arg <- function(eta,beta,lambda,y,sigma2,alpha)     # value of MY-envelope
{
  MY_env <- alpha*(sum(abs(D_mat%*%eta))) + sum((beta-eta)^2)/(2*lambda) + 
                              sum((y-eta)^2)/(2*sigma2)    
  return(MY_env)
}

# function calculates the value of the proximal function
prox_func <- function(beta,lambda,alpha,sigma2,k,grid)
{
  lambda_star <- lambda/sigma2
  betaval <- (beta+(lambda_star*y))/sqrt(1+lambda_star)
  lambdaval <- (lambda*alpha)/sqrt(1+lambda_star)
  temp = trendfilter(grid,betaval, k=k,lambda = lambdaval)$beta
  out <- (temp)/sqrt(1+lambda_star)
  return(as.vector(out))
}

log_gradpi <- function(beta,lambda,y,sigma2,alpha,k,grid)  # gradient of log target
{
  beta_prox <- prox_func(beta,lambda,alpha,sigma2,k,grid)
  ans <-  (beta-beta_prox)/lambda
  return(-ans)
}

# bark.prop <- function(val, delta, lambda)
# {
#   y <- rnorm(1, 0, sqrt(delta))
#   prob <- 1 / (1 + exp( -y*log_gradpi(val, lambda)))
#   ifelse(runif(1) <= prob, prop <- val + y, prop <- val - y)
#   return(prop)
# }
# 
# bark.dens <- function(in_val, propval, delta, lambda)
# {
#   numer <- 2 * dnorm(propval, in_val, sqrt(delta))
#   denom <- 1 + exp( - (propval - in_val) * log_gradpi(in_val, lambda))
#   value <- numer / denom
#   return(value)
# }

# MYMALA sampling for covariance matrix estimation

mymala_cov_fn <- function(y, alpha, sigma2, k, grid, iter, delta)
{
  samp.mym <- matrix(0, nrow = iter, ncol = length(y))
  lambda <- lamb_coeff
  beta_current <- trendfilter(grid,y, k=k,lambda = sigma2*alpha)$beta
  samp.mym[1,] <- beta_current
  accept <- 0
  for (i in 2:iter) 
  {
    beta_next <- rnorm(length(beta_current), beta_current + 
                         (delta / 2)*log_gradpi(beta_current,lambda,y,sigma2,alpha,k,grid), 
                       sqrt(delta))   # proposal step
    prox_val.next <- prox_func(beta_next, lambda, alpha, sigma2, k, grid)
    prox_val.curr <- prox_func(beta_current, lambda, alpha,sigma2, k, grid)
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
      accept <- accept + 1
    }
    else
    {
      samp.mym[i,] <- beta_current
    }
    beta_current <- samp.mym[i,]
    if(i %% 1000 == 0){
      print(i)
    }
  }
  print(accept/iter)
  object <- samp.mym
  cov_mat <- cov(object)
  result <- list(object, cov_mat)
  return(result)
}

dmvnorm_fn <- function(point, mu, mat, delta)
{
  diff <- point - mu
  exp_term <- (t(diff) %*% mat) %*% diff 
  den_value <- -(exp_term/(2*delta))
}

mymala <- function(y, alpha, sigma2, k, grid, iter, delta, covmat)
{
  samp.mym <- matrix(0, nrow = iter, ncol = length(y))
  lambda <- lamb_coeff
  wts_is_est <- numeric(length = iter)
  beta_current <- trendfilter(grid,y, k=k,lambda = sigma2*alpha)$beta
  samp.mym[1,] <- beta_current
  g_val <- alpha*sum(abs(D_mat%*%beta_current)) + sum((y - beta_current)^2)/(2*sigma2)
  prox_start <- prox_func(beta_current, lambda, alpha, sigma2, k, grid)
  g_lambda_val <- prox_arg(prox_start, beta_current, lambda=lambda, y, sigma2, alpha)
  wts_is_est[1] <- g_lambda_val - g_val
  U <- sqrtm(covmat)
  mat.inv <- solve(covmat)
  accept <- 0
  for (i in 2:iter) 
  {
    beta_next <- beta_current +  ((delta / 2)*covmat)%*%log_gradpi(beta_current,lambda,y,sigma2,alpha,k,grid) + 
      (sqrt(delta)*U) %*% rnorm(length(beta_current), 0, 1)   # proposal step
    prox_val.next <- prox_func(beta_next, lambda, alpha, sigma2, k, grid)
    prox_val.curr <- prox_func(beta_current, lambda, alpha, sigma2, k, grid)
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
    if(i %% 1000 == 0){
      print(i)
    }
  }
  print(accept/iter)
  object <- list(samp.mym, wts_is_est)
  return(object)
}

px.mala <- function(y, alpha, sigma2, k, grid, iter, delta, covmat)
{
  samp.pxm <- matrix(0, nrow = iter, ncol = length(y))
  lambda <- lamb_coeff*sigma2
  beta_current <- trendfilter(grid,y, k=k,lambda = sigma2*alpha)$beta
  samp.pxm[1,] <- beta_current
  accept <- 0
  U <- sqrtm(covmat)
  mat.inv <- solve(covmat)
  for (i in 2:iter) 
  {
    beta_next <- beta_current +  ((delta / 2)*covmat)%*%log_gradpi(beta_current,lambda,y,sigma2,alpha,k,grid) + 
      (sqrt(delta)*U) %*% rnorm(length(beta_current), 0, 1)   # proposal step
    U_betanext <- - (sum((y - beta_next)^2)/(2*sigma2) + alpha*(sum(abs(D_mat%*%beta_next))))
    U_betacurr <- - (sum((y - beta_current)^2)/(2*sigma2) + alpha*(sum(abs(D_mat%*%beta_current))))
    q.next_to_curr <- dmvnorm_fn(beta_current, beta_next + 
                                   ((delta / 2)*covmat)%*%log_gradpi(beta_next,lambda,y,sigma2,alpha,k,grid), mat.inv, delta)
    q.curr_to_next <- dmvnorm_fn(beta_next, beta_current + 
                                   ((delta / 2)*covmat)%*%log_gradpi(beta_current,lambda,y,sigma2,alpha,k,grid), mat.inv, delta) 
    mh.ratio <- U_betanext + q.next_to_curr - (U_betacurr + q.curr_to_next)
    if(log(runif(1)) <= mh.ratio)
    {
      samp.pxm[i,] <- beta_next
      accept <- accept + 1
    }
    else
    {
      samp.pxm[i,] <- beta_current
    }
    beta_current <- samp.pxm[i,]
    if(i %% 1000 == 0){
      print(i)
    }
  }
  print(accept/iter)
  return(samp.pxm)
}


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
