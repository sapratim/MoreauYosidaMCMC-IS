## Code for asymptotic variance comparison of weighted importance sampling estimator 
## for trend filtering example
library(mcmcse)
library(coda)
library(glmgen)
library(Matrix)
library(spmrf)

alpha_hat <- 18.9   # obtained from the first dataset
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

D_mat <- getD(k=1, n=1e2, x)   #  D matrix

log_target <- function(eta,beta,lambda,y,sigma2,alpha)
{
  f.beta <- sum((y-beta)^2)/(2*sigma2)
  g_lambda.beta <- alpha*(sum(abs(D_mat%*%eta))) + sum((beta-eta)^2)/(2*lambda)
  dens_val <- f.beta + g_lambda.beta
  return(-dens_val)
}

prox_arg <- function(eta,beta,lambda,alpha)     # MY-envelope
{
  MY_env <- alpha*(sum(abs(D_mat%*%eta))) + sum((beta-eta)^2)/(2*lambda)    
  return(MY_env)
}

# function calculates the value of the proximal function
prox_func <- function(beta,lambda,alpha,k,grid)
{
  out = trendfilter(grid,beta, k=k,lambda = lambda*alpha)
  return(as.vector(out$beta))
}

log_gradpi <- function(beta,lambda,y,sigma2,alpha,k,grid)
{
  beta_prox <- prox_func(beta,lambda,alpha,k,grid)
  ans <- (beta-y)/sigma2 + (beta-beta_prox)/lambda
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

# MYMALA sampling and importance sampling weights

mymala <- function(y, alpha, sigma2, k, grid, iter, delta)
{
  samp.mym <- matrix(0, nrow = iter, ncol = length(y))
  lambda <- 0.00001*sigma2
  wts_is_est <- numeric(length = iter)
  beta_current <- y
  samp.mym[1,] <- beta_current
  g_val <- alpha*sum(abs(D_mat%*%beta_current))
  prox_start <- prox_func(beta_current, lambda, alpha, k, grid)
  g_lambda_val <- prox_arg(prox_start, beta_current, lambda=lambda, alpha)
  wts_is_est[1] <- exp(g_lambda_val - g_val)
  accept <- 0
  for (i in 2:iter) 
  {
    beta_next <- rnorm(length(beta_current), beta_current + 
            (delta / 2)*log_gradpi(beta_current,lambda,y,sigma2,alpha,k,grid), 
            sqrt(delta))   # proposal step
    prox_val.next <- prox_func(beta_next, lambda, alpha, k, grid)
    prox_val.curr <- prox_func(beta_current, lambda, alpha, k, grid)
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
      g_val <- alpha*sum(abs(D_mat%*%beta_next))
      g_lambda_val <- prox_arg(prox_val.next, beta_next, lambda=lambda, alpha)
      wts_is_est[i] <- exp(g_lambda_val - g_val)
      accept <- accept + 1
    }
    else
    {
      samp.mym[i,] <- beta_current
      g_val <- alpha*sum(abs(D_mat%*%beta_current))
      g_lambda_val <- prox_arg(prox_val.curr, beta_current, lambda=lambda, alpha)
      wts_is_est[i] <- exp(g_lambda_val - g_val)
    }
    beta_current <- samp.mym[i,]
  }
  object <- list(samp.mym, wts_is_est)
  return(object)
}

px.mala <- function(y, alpha, sigma2, k, grid, iter, delta)
{
  samp.pxm <- matrix(0, nrow = iter, ncol = length(y))
  lambda <- 0.00001*sigma2
  beta_current <- y
  samp.pxm[1,] <- beta_current
  accept <- 0
  for (i in 2:iter) 
  {
    beta_next <- rnorm(length(beta_current), beta_current + 
                         (delta / 2)*log_gradpi(beta_current,lambda,y,sigma2,alpha,k,grid), 
                       sqrt(delta))   # proposal step
    U_betanext <- - (sum((y - beta_next)^2)/(2*sigma2) + alpha*(sum(abs(D_mat%*%beta_next))))
    U_betacurr <- - (sum((y - beta_current)^2)/(2*sigma2) + alpha*(sum(abs(D_mat%*%beta_current))))
    q.next_to_curr <- sum(dnorm(beta_current, beta_next + 
                                  (delta / 2)*log_gradpi(beta_next,lambda,y,sigma2,alpha,k,grid),
                                sqrt(delta), log = TRUE))
    q.curr_to_next <- sum(dnorm(beta_next, beta_current + 
                                  (delta / 2)*log_gradpi(beta_current,lambda,y,sigma2,alpha,k,grid),
                                sqrt(delta), log = TRUE))
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
  }
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

iter <- 1e3
delta <- 0.0008
mymala(y, alpha_hat, sigma2_hat, k=1, grid=x, iter, delta)
px.mala(y, alpha_hat, sigma2_hat, k=1, grid=x, iter, delta)


 asymp_covmat_is <- matrix(0, length(y), length(y))
 asymp_covmat_pxm <- matrix(0, length(y), length(y))
 
  mala.is <- mymala(y, alpha_hat, sigma2_hat, k=1, grid=x, iter, delta)
  px_mala <- px.mala(y, alpha_hat, sigma2_hat, k=1, grid=x, iter, delta)

   # Importance sampling estimator asymptotic variance

   is_samp <- as.matrix(unlist(mala.is[1]))
   is_wts <- as.numeric(unlist(mala.is[2]))
   wts_mean <- mean(is_wts)
   num <- is_samp*is_wts
   sum_mat <- apply(num, 2, sum)
   is_est <- sum_mat / sum(is_wts)
   input_mat <- cbind(num, is_wts)  # input samples for mcse
   Sigma_mat <- mcse.multi(input_mat)$cov  # estimated covarince matrix of the tuple
   kappa_eta_mat <- cbind(diag(1/wts_mean, length(y)), t(is_est/wts_mean)) # derivative of kappa matrix
   asymp_covmat_is <- (kappa_eta_mat %*% Sigma_mat) %*% t(kappa_eta_mat)

  # PxMALA asymptotic variance

   asymp_covmat_pxm[i] <- mcse.multi(px_mala)$cov

   # PxBarker asymptotic variance

   # asymp_covmat_pxb[i] <- mcse.multi(px_bark)$cov


 # Asymptotic variance comparison
 var_mat <- cbind(det(asymp_covmat_is), det(asymp_covmat_pxm)) #asymp_covmat_pxb)
 colnames(var_mat) <- c("Imp_sampling", "PxMala")

 var_mat
