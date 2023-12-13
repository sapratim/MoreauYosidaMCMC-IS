
## Code for asymptotic variance comparison of weighted importance sampling estimator 
## for trend filtering example
library(mcmcse)
library(coda)
library(glmgen)
library(Matrix)
library(spmrf)
library(expm)
load("MC_pcm.Rdata")

set.seed(12345)
alpha_hat <- 5   # obtained from the first dataset
sigma2_hat <- 1  # obtained from the first dataset
x <- seq(1,100,len=100)
 f <- Vectorize(function(x){if(x<=35){x} else if(x<=70){70-x} else{0.5*x-35}})
fx_linear <- f(x)
y <- fx_linear + rnorm(length(x), sd = 1)
tol <- 1e-20
tol_pcm <- 1e-8
max_iter <- 1e5L
max_iter_pcm <- 1e3L

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
  betaval <- (beta*sigma2+(lambda*y))/(sigma2+lambda)
  lambdaval <- (alpha)/ ((lambda + sigma2)/ (lambda*sigma2) )
  temp = trendfilter(grid,betaval, k=k,lambda = lambdaval,
                control = trendfilter.control.list(obj_tol = tol, max_iter = max_iter))$beta
  return(as.vector(temp))
}

prox_func_pcm <- function(beta,lambda,alpha,sigma2,k,grid)
{
  betaval <- (beta*sigma2+(lambda*y))/(sigma2+lambda)
  lambdaval <- (alpha)/ ((lambda + sigma2)/ (lambda*sigma2) )
  temp = trendfilter(grid,betaval, k=k,lambda = lambdaval,
           control = trendfilter.control.list(obj_tol = tol_pcm, max_iter = max_iter_pcm))$beta
  return(as.vector(temp))
}

log_gradpi <- function(beta,lambda,y,sigma2,alpha,k,grid)  # gradient of log target
{
  beta_prox <- prox_func(beta,lambda,alpha,sigma2,k,grid)
  ans <-  (beta-beta_prox)/lambda
  return(-ans)
}

dmvnorm_fn <- function(point, mu, mat, delta)
{
  diff <- point - mu
  exp_term <- (t(diff) %*% mat) %*% diff 
  den_value <- -(exp_term/(2*delta))
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

