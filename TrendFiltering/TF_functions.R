############################################
## Main functions for trendfiltering example
## Plus data generation
############################################

# loading libraries needed
library(mcmcse)
library(glmgen)
library(Matrix)
library(expm)
library(foreach)
library(doParallel)

# generating data
set.seed(12345)
alpha_hat <- 5   # obtained from the first dataset
sigma2_hat <- 9  # obtained from the first dataset
x <- seq(1,100,len=100)
f <- Vectorize(function(x){if(x<=35){x} else if(x<=70){70-x} else{0.5*x-35}})
fx_linear <- f(x)
y <- fx_linear + rnorm(length(x), sd = 3)


# setting parameters for the proximal mappings
tol <- 1e-6
max_iter <- 200L


# function calculates the D matrix of the penalty term
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

# log target of pi-lambda
log_pilambda <- function(eta,beta,lambda,y,sigma2,alpha)
{
  dens_val <- alpha*(sum(abs(D_mat%*%eta))) + sum((beta-eta)^2)/(2*lambda) + 
    sum((y-eta)^2)/(2*sigma2)
  return(-dens_val)
}

log_pi <- function(beta,y,sigma2,alpha)
{
  dens_val <- alpha*(sum(abs(D_mat%*%beta))) + sum((y-beta)^2)/(2*sigma2)
  return(-dens_val)
}


# value of MY-envelope
prox_arg <- function(eta,beta,lambda,y,sigma2,alpha)     
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


# gradient of log target (pi-lambda)
grad_logpiLam <- function(beta,lambda,y,sigma2,alpha,k,grid)  
{
  beta_prox <- prox_func(beta,lambda,alpha,sigma2,k,grid)
  ans <-  (beta-beta_prox)/lambda
  return(-ans)
}

# evaluate dnorm for multivariate normal
dmvnorm_fn <- function(point, mu, mat, delta)
{
  diff <- point - mu
  exp_term <- (t(diff) %*% mat) %*% diff 
  den_value <- -(exp_term/(2*delta))
}

bark.prop <- function(beta,lambda,y,sigma2,alpha,k,grid,delta)
{
  aux_var <- rnorm(length(beta), 0, 1)
  y <- delta*aux_var
  denom_prod <- y*grad_logpiLam(beta,lambda,y,sigma2,alpha,k,grid)
  prob <- 1 / (1 + exp(- sum(denom_prod)))
  ifelse(runif(1) <= prob, prop <- beta + y, prop <- beta - y)
  return(prop)
}

log_q_ratio_barker<-function(x,y,grad_x,grad_y)
{
  # x: current location (vector)
  # y: proposed location (vector)
  # grad_x: target log-posterior gradient at x (vector)
  # grad_y: target log-posterior gradient at y (vector)
  a1<-  c(-(x-y)*grad_y)
  a2<-  c(-grad_x*(y-x))
  return(sum(-(pmax(a1,0)+log1p(exp(-abs(a1))))+
               (pmax(a2,0)+log1p(exp(-abs(a2))))
  ))
}

quant <- function(j, mat)     ### function for quantile estimation
{
  mat_ordered <- mat[order(mat[,j], decreasing = FALSE), ]
  order_comp <- mat_ordered[,j]
  weights_order <- mat_ordered[, length(y)+1]
  bound_mat <- cbind(order_comp, weights_order)
  return(bound_mat)
}


#### MALA for pi-lambda without weights, for warm-up

mymala_cov_fn <- function(y, alpha, sigma2, k, grid, iter, delta) # pre-conditioned mala
{
  nvar <- length(y)
  samp.mym <- matrix(0, nrow = iter, ncol = nvar)
  lambda <- lamb_coeff
  
  # starting value computations
  beta_current <- y
  prox_val.curr <- prox_func(beta_current, lambda, alpha,sigma2, k, grid)
  targ_val.curr <- log_pilambda(prox_val.curr,beta_current,lambda,y,sigma2,alpha)
  samp.mym[1,] <- beta_current
  accept <- 0
  
  # mcmc
  for (i in 2:iter) 
  {
    
    # proposal
    prop.mean <- beta_current + (delta / 2)*grad_logpiLam(beta_current,lambda,y,sigma2,alpha,k,grid)
    beta_next <- rnorm(nvar, prop.mean, sd = sqrt(delta))
    
    # Evaluating proximal map and target (pi-lambda)
    prox_val.next <- prox_func(beta_next, lambda, alpha, sigma2, k, grid)
    targ_val.next <- log_pilambda(prox_val.next,beta_next,lambda,y,sigma2,alpha)
    
    other.mean <-  beta_next + (delta / 2)*grad_logpiLam(beta_next,lambda,y,sigma2,alpha,k,grid)
    q.next_to_curr <- sum(dnorm(beta_current, other.mean, sd = sqrt(delta), log = TRUE))
    q.curr_to_next <- sum(dnorm(beta_next, prop.mean, sd = sqrt(delta), log = TRUE))
    
    mh.ratio <- targ_val.next + q.next_to_curr - (targ_val.curr + q.curr_to_next)  # mh  ratio
    # print(mh.ratio)
    if(log(runif(1)) <= mh.ratio)
    {
      samp.mym[i,] <- beta_next
      prox_val.curr <- prox_val.next
      targ_val.curr <- targ_val.next
      accept <- accept + 1
    }
    else
    {
      samp.mym[i,] <- beta_current
    }
    beta_current <- samp.mym[i,]
    if(i %% 1000 == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(accept/iter)
  object <- samp.mym
  return(object)
}

##### MYMALA samples function

mymala <- function(y, alpha, sigma2, k, grid, iter, delta, start)
{
  nvar <- length(y)
  samp.mym <- matrix(0, nrow = iter, ncol = nvar)
  lambda <- lamb_coeff
  wts_is_est <- numeric(length = iter)
  
  # starting value computations
  beta_current <- start
  prox_val.curr <- prox_func(beta_current, lambda, alpha,sigma2, k, grid)
  targ_val.curr <- log_pilambda(prox_val.curr,beta_current,lambda,y,sigma2,alpha)
  samp.mym[1,] <- beta_current
  accept <- 0
  
  # weights calculation
  g_val <- -log_pi(beta_current, y, sigma2,alpha) 
  g_lambda_val <- -targ_val.curr
  wts_is_est[1] <- g_lambda_val - g_val
  
  # for pre-conditioned MALA
  mat.inv <- diag(1, nvar, nvar)
  accept <- 0
  for (i in 2:iter) 
  {
    # proposal step
    prop.mean <- beta_current + (delta / 2)*grad_logpiLam(beta_current,lambda,y,sigma2,alpha,k,grid)
    beta_next <-  prop.mean + (sqrt(delta))*rnorm(nvar, 0, 1)  
    
    # calculating prox
    prox_val.next <- prox_func(beta_next, lambda, alpha, sigma2, k, grid)
    targ_val.next <- log_pilambda(prox_val.next,beta_next,lambda,y,sigma2,alpha)
    
    q.next_to_curr <- dmvnorm_fn(beta_current, beta_next + 
                                   (delta / 2)*grad_logpiLam(beta_next,lambda,y,sigma2,alpha,k,grid),
                                 mat.inv, delta)
    q.curr_to_next <- dmvnorm_fn(beta_next, prop.mean, mat.inv, delta) 
    
    
    mh.ratio <- targ_val.next + q.next_to_curr - (targ_val.curr + q.curr_to_next)  # mh  ratio
    # print(mh.ratio)
    if(log(runif(1)) <= mh.ratio)
    {
      samp.mym[i,] <- beta_next
      targ_val.curr <- targ_val.next
      prox_val.curr <- prox_val.next
      
      # weights
      g_val <- -log_pi(beta_next, y, sigma2,alpha)
      g_lambda_val <- - targ_val.curr
      wts_is_est[i] <- g_lambda_val - g_val
      accept <- accept + 1
    }
    else
    {
      samp.mym[i,] <- beta_current
      g_val <- -log_pi(beta_current, y, sigma2,alpha)
      g_lambda_val <- - targ_val.curr
      wts_is_est[i] <- g_lambda_val - g_val
    }
    beta_current <- samp.mym[i,]
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(accept/iter)
  object <- list(samp.mym, wts_is_est)
  return(object)
}

##### PxMALA samples function

px.mala <- function(y, alpha, sigma2, k, grid, iter, delta, start)
{
  nvar <- length(y)
  samp.pxm <- matrix(0, nrow = iter, ncol = nvar)
  lambda <- lamb_coeff
  
  # starting value computations
  beta_current <- start
  samp.pxm[1,] <- beta_current
  U_betacurr <- log_pi(beta_current, y, sigma2, alpha)
  accept <- 0
  
  mat.inv <- diag(1, nvar, nvar)
  for (i in 2:iter)
  {
    # proposal step
    prop.mean <- beta_current +  (delta / 
                                    2)*grad_logpiLam(beta_current,lambda,y,sigma2,alpha,k,grid)
    beta_next <-  prop.mean + sqrt(delta)*rnorm(nvar, 0, 1) 
    
    U_betanext <- log_pi(beta_next, y, sigma2, alpha)
    q.next_to_curr <- dmvnorm_fn(beta_current, beta_next +
                                   (delta /  2)*grad_logpiLam
                                 (beta_next,lambda,y,sigma2,alpha,k,grid), 
                                 mat.inv, delta)
    q.curr_to_next <- dmvnorm_fn(beta_next, prop.mean, mat.inv, delta)
    
    mh.ratio <- U_betanext + q.next_to_curr - (U_betacurr + q.curr_to_next) # mh ratio
    
    if(log(runif(1)) <= mh.ratio)
    {
      samp.pxm[i,] <- beta_next
      U_betacurr <- U_betanext
      accept <- accept + 1
    }
    else
    {
      samp.pxm[i,] <- beta_current
    }
    beta_current <- samp.pxm[i,]
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(accept/iter)
  return(samp.pxm)
}

## MYBarker samples function

mybarker <- function(y, alpha, sigma2, k, grid, iter, delta, start)
{
  nvar <- length(y)
  samp.bark <- matrix(0, nrow = iter, ncol = nvar)
  lambda <- lamb_coeff
  wts_is_est <- numeric(length = iter)
  
  # starting value computations
  beta_current <- start
  samp.bark[1,] <- beta_current
  prox_val.curr <- prox_func(beta_current, lambda, alpha, sigma2, k, grid)
  targ_val.curr <- log_pilambda(prox_val.curr,beta_current,lambda,y,sigma2,alpha)
  accept <- 0
  
  # weights calculation
  g_val <- - log_pi(beta_current, y, sigma2, alpha)
  g_lambda_val <- - targ_val.curr
  wts_is_est[1] <- g_lambda_val - g_val
  
  # using pre conditioned mala covariance
  mat.inv <- diag(1, nvar, nvar)
  
  for (i in 2:iter) 
  {
    # proposal step
    beta_next <- bark.prop(beta_current,lambda,y,sigma2,alpha,k,grid,delta)
    
    prox_val.next <- prox_func(beta_next, lambda, alpha, sigma2, k, grid)
    targ_val.next <- log_pilambda(prox_val.next,beta_next,lambda,y,sigma2,alpha)
    
    grad_beta_curr <- grad_logpiLam(beta_current, lambda, y, sigma2, alpha, k, grid)
    grad_beta_next <- grad_logpiLam(beta_next, lambda, y, sigma2, alpha, k, grid)
    
    mh.ratio <- targ_val.next - targ_val.curr + 
      log_q_ratio_barker(beta_current,beta_next,grad_beta_curr,grad_beta_next)
    
    if(log(runif(1)) <= mh.ratio)
    {
      samp.bark[i,] <- beta_next
      prox_val.curr <- prox_val.next
      targ_val.curr <- targ_val.next
      
      # weights
      g_val <- - log_pi(beta_next, y, sigma2,alpha)
      g_lambda_val <- - targ_val.curr
      wts_is_est[i] <- g_lambda_val - g_val
      accept <- accept + 1
    }
    else
    {
      samp.bark[i,] <- beta_current
      g_val <- - log_pi(beta_current, y, sigma2,alpha)
      g_lambda_val <- - targ_val.curr
      wts_is_est[i] <- g_lambda_val - g_val
    }
    beta_current <- samp.bark[i,]
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(acc_rate <- accept/iter)
  object <- list(samp.bark, wts_is_est, acc_rate)
  return(object)
}

## PxBarker samples

px.barker <- function(y, alpha, sigma2, k, grid, iter, delta, start)
{
  nvar <- length(y)
  samp.bark <- matrix(0, nrow = iter, ncol = length(y))
  lambda <- lamb_coeff
  
  # starting value computations
  beta_current <- start
  samp.bark[1,] <- beta_current
  U_betacurr <- log_pi(beta_current, y, sigma2, alpha)
  
  # using pre conditioned mala covariance
  mat.inv <- diag(1, nvar, nvar)
  accept <- 0
  for (i in 2:iter)
  {
    # proposal step
    beta_next <- bark.prop(beta_current,lambda,y,sigma2,alpha,k,grid,delta)
    
    U_betanext <- log_pi(beta_next, y, sigma2, alpha)
    grad_beta_curr <- grad_logpiLam(beta_current, lambda, y, sigma2, alpha, k, grid)
    grad_beta_next <- grad_logpiLam(beta_next, lambda, y, sigma2, alpha, k, grid)
    
    mh.ratio <- U_betanext - U_betacurr + 
      log_q_ratio_barker(beta_current,beta_next,grad_beta_curr,grad_beta_next)
    if(log(runif(1)) <= mh.ratio)
    {
      samp.bark[i,] <- beta_next
      U_betacurr <- U_betanext
      accept <- accept + 1
    }
    else
    {
      samp.bark[i,] <- beta_current
    }
    beta_current <- samp.bark[i,]
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(acc_rate <- accept/iter)
  object <- list(samp.bark, acc_rate)
  return(object)
}

##  myhmc samples

myhmc <- function(y, alpha, sigma2, k, grid, iter, eps_hmc, L, start)
{
  nvar <- length(y)
  samp.hmc <- matrix(0, nrow = iter, ncol = nvar)
  lambda <- lamb_coeff
  wts_is_est <- numeric(length = iter)
  
  # starting value computations
  beta <- start
  proxval_curr <- prox_func(beta, lambda, alpha, sigma2, k, grid)
  samp.hmc[1,] <- beta
  
  # weights calculation
  g_val <- - log_pi(beta, y, sigma2, alpha)
  g_lambda_val <- - log_pilambda(proxval_curr, beta, lambda=lambda, y, sigma2, alpha)
  wts_is_est[1] <- g_lambda_val - g_val
  
  # proposal step
  mom_mat <- matrix(rnorm(iter*length(y)), nrow = iter, ncol = length(y))
  accept <- 0
  for (i in 2:iter) 
  {
    p_prop <- mom_mat[i,]
    U_beta <- -grad_logpiLam(beta, lambda,y,sigma2,alpha,k,grid)
    p_current <- p_prop - eps_hmc*U_beta /2  # half step for momentum
    q_current <- beta
    for (j in 1:L)
    {
      beta <- beta + eps_hmc*p_current   # full step for position
      U_beta <- -grad_logpiLam(beta, lambda,y,sigma2,alpha,k,grid)
      if(j!=L) p_current <- p_current - eps_hmc*U_beta  # full step for momentum
    }
    p_current <- p_current - eps_hmc*U_beta/2
    p_current <- - p_current  # negation to make proposal symmetric
    
    #  calculating prox
    proxval_prop <- prox_func(beta = beta,lambda,alpha,sigma2,k,grid)
    U_curr <- - log_pilambda(proxval_curr,q_current,lambda,y,sigma2,alpha)
    U_prop <- - log_pilambda(proxval_prop,beta,lambda,y,sigma2,alpha)
    K_curr <-  sum((p_prop^2)/2)
    K_prop <-  sum((p_current^2)/2)
    
    log_acc_prob = U_curr - U_prop + K_curr - K_prop    # mh ratio
    
    if(log(runif(1)) <= log_acc_prob )
    {
      samp.hmc[i,] <- beta
      proxval_curr <- proxval_prop
      
      # weights
      g_val <- -log_pi(beta, y, sigma2, alpha)
      g_lambda_val <- - log_pilambda(proxval_prop, beta, lambda=lambda, y, sigma2, alpha)
      wts_is_est[i] <- g_lambda_val - g_val
      accept <- accept + 1
    }
    else
    {
      samp.hmc[i,] <- q_current
      g_val <- -log_pi(q_current, y, sigma2, alpha)
      g_lambda_val <- - log_pilambda(proxval_curr, q_current, lambda=lambda, y, sigma2, alpha)
      wts_is_est[i] <- g_lambda_val - g_val
      beta <- q_current
    }
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))}
  } 
  print(acc_rate <- accept/iter)
  object <- list(samp.hmc, wts_is_est, acc_rate)
  return(object)
}

## pxhmc samples

pxhmc <- function(y, alpha, sigma2, k, grid, iter, eps_hmc, L, start)
{
  nvar <- length(y)
  samp.hmc <- matrix(0, nrow = iter, ncol = nvar)
  lambda <- lamb_coeff
  
  # starting value computations
  beta <- start
  samp.hmc[1,] <- beta
  
  # proposal step
  mom_mat <- matrix(rnorm(iter*length(y)), nrow = iter, ncol = length(y))
  accept <- 0
  
  for (i in 2:iter) 
  {
    p_prop <- mom_mat[i,]
    U_beta <- -grad_logpiLam(beta, lambda,y,sigma2,alpha,k,grid)
    p_current <- p_prop - eps_hmc*U_beta /2  # half step for momentum
    q_current <- beta
    for (j in 1:L)
    {
      beta <- beta + eps_hmc*p_current   # full step for position
      U_beta <- -grad_logpiLam(beta, lambda,y,sigma2,alpha,k,grid)
      if(j!=L) p_current <- p_current - eps_hmc*U_beta  # full step for momentum
    }
    p_current <- p_current - eps_hmc*U_beta/2
    p_current <- - p_current  # negation to make proposal symmetric
    
    U_curr <- - log_pi(q_current, y, sigma2, alpha)
    U_prop <- - log_pi(beta, y, sigma2, alpha)
    K_curr <-  sum((p_prop^2)/2)
    K_prop <-  sum((p_current^2)/2)
    
    log_acc_prob = U_curr - U_prop + K_curr - K_prop
    
    if(log(runif(1)) <= log_acc_prob )
    {
      samp.hmc[i,] <- beta
      accept <- accept + 1
    }
    else
    {
      samp.hmc[i,] <- q_current
      beta <- q_current
    }
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))}
  } 
  print(acc_rate <- accept/iter)
  object <- list(samp.hmc, acc_rate)
  return(object)
}
