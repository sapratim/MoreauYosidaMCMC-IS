# R Code for Proximal of the Nuclear-Norm problem:
# softthreshold(seq(-10, 10, length.out=10), 1) 

library(mcmcse)
library(Matrix)
library(expm)
library(foreach)
library(doParallel)
library(ks)

set.seed(8024248)
n <- 64
a <- 8
mat <- matrix(0, nrow = n, ncol = n)
vec.mat <- rep(c(1, 0), each = a, times = n/(2*a))
checker <- matrix(0, nrow = n, ncol = n)
for(j in seq(1, n/2, by = a))
{
  for(k in 1:a)
  {
    checker[j+k-1, ] <- vec.mat
  }
  vec.mat <- rev(vec.mat)
}
for(j in seq((n/2+1), n, by = a))
{
  for(k in 1:a)
  {
    checker[j+k-1, ] <- (vec.mat == 0)*(vec.mat)  + (vec.mat == 1)*(vec.mat - .50)
  }
  vec.mat <- rev(vec.mat)
}
noise <- matrix(rnorm(n^2, 0, sqrt(0.01)), nrow = n, ncol = n)
image_mat <- checker + noise

x <- vec(checker)
y <- vec(image_mat)
######### Display the checkerboard matrix as an image

# image(checker, col = gray.colors(4, start = 0, end = 1), axes = FALSE)
# image(image_mat, col = gray.colors(4, start = 0, end = 1), axes = FALSE)

##############------Functions------##################################

nucl_norm <- function(vect)    ## vector input
{
  A <- matrix(vect, nrow = n, ncol = n)
  norm_val <- sum(svd(A)$d)
  return(norm_val)
}

log_pi <- function(x,y,sigma2,alpha)
{
  n_norm <- nucl_norm(x)
  dens_val <- alpha*n_norm + sum((y - x)^2)/(2*sigma2)
  return(-dens_val)
}

log_pilambda <- function(eta,x,lambda,y,sigma2,alpha)   # log target function of pi^lambda
{
  n_norm <- nucl_norm(eta)
  dens_val <- alpha*n_norm + sum((eta-x)^2)/(2*lambda) + 
    sum((y-eta)^2)/(2*sigma2)
  return(-dens_val)
}

#########  Soft threshold function

softthreshold <- function(u, lambda) {       ####  u is a vector
  return(sign(u)*sapply(u, FUN=function(x) {max(abs(x)-lambda,0)}))
}

######### Proximity mapping

prox_func <- function(x,lambda,y,sigma2,alpha) {   #### input x and y as a vector
  num <- lambda*y + sigma2*x
  denom <- lambda + sigma2
  mat <- matrix(num/denom, nrow = n, ncol = n)
  svdsol <- svd(mat)
  s <- softthreshold(svdsol$d, (alpha*sigma2*lambda)/denom)
  output <- svdsol$u %*% (s*t(svdsol$v))  # Multiply each row of V’ by singular values
  return(vec(output))
}

grad_logpiLam <- function(x,lambda,y,sigma2,alpha)  # gradient of log target
{
  x_prox <- prox_func(x,lambda,y,sigma2,alpha)
  ans <-  (x-x_prox)/lambda
  return(-ans)
}

# function for barker proposal
bark.prop <- function(x, alpha, lambda, y, sigma2, delta)
{
  aux_var <- rnorm(length(x), 0, 1)
  z <- sqrt(delta)*aux_var
  denom_prod <- z*grad_logpiLam(x,lambda,y,sigma2,alpha)
  prob <- 1 / (1 + exp(- denom_prod))
  unifs <- runif(length(x))
  prop <- (x + z)*(unifs <= prob) + (x - z)*(unifs > prob)
  return(prop)
}

# barker log density
log_bark.dens <- function(curr_point, prop_point, grad_curr_point, delta)
{
  rw_dens <- sum(dnorm(prop_point-curr_point, 0, sqrt(delta), log = TRUE))
  exp_term <- - (grad_curr_point*(prop_point-curr_point))
  denom <- sum(log1p(exp(exp_term)))
  dens_val <- rw_dens - denom
  return(dens_val)
}

# Asymptotic covariance matrix

asymp_cov_func <- function(chain, weights)
{
  dim <- ncol(chain)
  wts_mean <- mean(weights)
  asymp_var<- numeric(length = dim)
  for (j in 1:dim) 
  {
    num <- chain[, j]*weights
    sum_mat <- sum(num)
    is_est <- sum_mat / sum(weights)
    input_mat <- cbind(num, exp(weights))  # input samples for mcse
    Sigma_mat <- mcse.multi(input_mat)$cov  # estimated covariance matrix of the tuple
    kappa_eta_mat <- cbind(1/wts_mean, -is_est/wts_mean) # derivative of kappa matrix
    asymp_var[j] <- (kappa_eta_mat %*% Sigma_mat) %*% t(kappa_eta_mat) # IS asymptotic variance
  }
  return(asymp_var)
}

##### Posterior mean
post_mean_fn <- function(chain, weights)
{
  chain_length <- nrow(chain)
  weight_mat <- matrix(0, nrow = chain_length, ncol = ncol(chain))
  for (i in 1:chain_length) {
    weight_mat[i,] <- chain[i,]*exp(weights[i])
  }
  num_sum <- apply(weight_mat, 2, sum)
  weights_sum <- sum(exp(weights))
  post_mean <- num_sum/weights_sum
  return(post_mean)
}


##### MYMALA samples function

mymala <- function(y, alpha, lambda, sigma2, iter, delta, start)
{
  nvar <- length(y)
  samp.mym <- matrix(0, nrow = iter, ncol = nvar)
  wts_is_est <- numeric(length = iter)
  
  # starting value computations
  samp_current <- start
  samp.mym[1,] <- samp_current
  prox_val.curr <- prox_func(samp_current, lambda, y, sigma2, alpha)
  targ_val.curr <- log_pilambda(prox_val.curr,samp_current,lambda,y,sigma2,alpha)
  
  # weights calculation
  psi_val <- - log_pi(samp_current,y,sigma2,alpha)
  psi_lambda_val <- - log_pilambda(prox_val.curr,samp_current,lambda,y,sigma2,alpha)
  wts_is_est[1] <- psi_lambda_val - psi_val
  
  accept <- 0
  for (i in 2:iter) 
  {
    # proposal step
    prop.mean <- samp_current + 
      (delta / 2)*grad_logpiLam(samp_current,lambda,y,sigma2,alpha)
    samp_next <- prop.mean + sqrt(delta)*rnorm(nvar, 0, 1)
    
    # calculating prox
    prox_val.next <- prox_func(samp_next,lambda,y,sigma2,alpha)
    targ_val.next <- log_pilambda(prox_val.next,samp_next,lambda,y,sigma2,alpha)
    
    q.next_to_curr <- sum(dnorm(samp_current, samp_next + 
                                  (delta / 2)*grad_logpiLam(samp_next,lambda,y,sigma2,alpha),
                                sqrt(delta), log = TRUE))
    q.curr_to_next <- sum(dnorm(samp_next, prop.mean,
                                sqrt(delta), log = TRUE))
    mh.ratio <- targ_val.next + q.next_to_curr - (targ_val.curr + q.curr_to_next)  # mh  ratio
    
    if(log(runif(1)) <= mh.ratio)
    {
      samp.mym[i,] <- samp_next
      prox_val.curr <- prox_val.next
      targ_val.curr <- targ_val.next
      
      # weights
      psi_val <- - log_pi(samp_next,y,sigma2,alpha)
      psi_lambda_val <- - targ_val.curr
      wts_is_est[i] <- psi_lambda_val - psi_val
      accept <- accept + 1
    }
    else
    {
      samp.mym[i,] <- samp_current
      psi_val <- - log_pi(samp_current,y,sigma2,alpha)
      psi_lambda_val <- - log_pilambda(prox_val.curr,samp_current,lambda,y,sigma2,alpha)
      wts_is_est[i] <- psi_lambda_val - psi_val
    }
    samp_current <- samp.mym[i,]
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

px.mala <- function(y, alpha, lambda, sigma2, iter, delta, start)
{
  nvar <- length(y)
  samp.pxm <- matrix(0, nrow = iter, ncol = nvar)
  
  # starting value computations
  samp_current <- start
  samp.pxm[1,] <- samp_current
  U_sampcurr <- log_pi(samp_current,y,sigma2,alpha)
  accept <- 0
  
  for (i in 2:iter)
  {
    # proposal step
    prop.mean <- samp_current + 
      (delta / 2)*grad_logpiLam(samp_current,lambda,y,sigma2,alpha)
    samp_next <- prop.mean + sqrt(delta)*rnorm(nvar, 0, 1)
    
    U_sampnext <- log_pi(samp_next,y,sigma2,alpha)
    q.next_to_curr <- sum(dnorm(samp_current, samp_next + 
                                  (delta / 2)*grad_logpiLam(samp_next,lambda,y,sigma2,alpha),
                                sqrt(delta), log = TRUE))
    q.curr_to_next <- sum(dnorm(samp_next, prop.mean,
                                sqrt(delta), log = TRUE))
    mh.ratio <- U_sampnext + q.next_to_curr - (U_sampcurr + q.curr_to_next)
    if(log(runif(1)) <= mh.ratio)
    {
      samp.pxm[i,] <- samp_next
      U_sampcurr <- U_sampnext
      accept <- accept + 1
    }
    else
    {
      samp.pxm[i,] <- samp_current
    }
    samp_current <- samp.pxm[i,]
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(accept/iter)
  return(samp.pxm)
}

## MYBarker samples function

mybarker <- function(y, alpha, lambda, sigma2, iter, delta, start)
{
  nvar <- length(y)
  samp.bark <- matrix(0, nrow = iter, ncol = nvar)
  wts_is_est <- numeric(length = iter)
  
  # starting value computations
  samp_current <- start
  samp.bark[1,] <- samp_current
  prox_val.curr <- prox_func(samp_current, lambda, y, sigma2, alpha)
  targ_val.curr <- log_pilambda(prox_val.curr,samp_current,lambda,y,sigma2,alpha)
  accept <- 0
  
  # weights calculation
  psi_val <- - log_pi(samp_current,y,sigma2,alpha)
  psi_lambda_val <- - log_pilambda(prox_val.curr,samp_current,lambda,y,sigma2,alpha)
  wts_is_est[1] <- psi_lambda_val - psi_val
  
  # For barker  
  for (i in 2:iter) 
  {
    # proposal step
    samp_next <- bark.prop(samp_current, alpha, lambda, y, sigma2, delta)
    
    prox_val.next <- prox_func(samp_next,lambda,y,sigma2,alpha)
    targ_val.next <- log_pilambda(prox_val.next,samp_next,lambda,y,sigma2,alpha)
    
    grad_samp_curr <- grad_logpiLam(samp_current,lambda,y,sigma2,alpha)
    grad_samp_next <- grad_logpiLam(samp_next,lambda,y,sigma2,alpha)
    
    mh.ratio <- targ_val.next + log_bark.dens(samp_next, samp_current, grad_samp_next, delta) - targ_val.curr - 
      log_bark.dens(samp_current, samp_next, grad_samp_curr, delta)
    if(log(runif(1)) <= mh.ratio)
    {
      samp.bark[i,] <- samp_next
      prox_val.curr <- prox_val.next
      targ_val.curr <- targ_val.next
      
      # weights
      psi_val <- - log_pi(samp_next, y, sigma2,alpha)
      psi_lambda_val <- - targ_val.curr
      wts_is_est[i] <- psi_lambda_val - psi_val
      accept <- accept + 1
    }
    else
    {
      samp.bark[i,] <- samp_current
      psi_val <- - log_pi(samp_current, y, sigma2,alpha)
      g_lambda_val <- - targ_val.curr
      wts_is_est[i] <- psi_lambda_val - psi_val
    }
    samp_current <- samp.bark[i,]
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

px.barker <- function(y, alpha, lambda, sigma2, iter, delta, start)
{
  nvar <- length(y)
  samp.bark <- matrix(0, nrow = iter, ncol = nvar)
  
  # starting value computations
  samp_current <- start
  samp.bark[1,] <- samp_current
  U_sampcurr <- log_pi(samp_current, y, sigma2, alpha)
  
  # For barker
  accept <- 0
  for (i in 2:iter)
  {
    # proposal step
    samp_next <- bark.prop(samp_current, alpha, lambda, y, sigma2, delta)
    
    U_sampnext <- log_pi(samp_next, y, sigma2, alpha)
    grad_samp_curr <- grad_logpiLam(samp_current,lambda,y,sigma2,alpha)
    grad_samp_next <- grad_logpiLam(samp_next,lambda,y,sigma2,alpha)
    
    mh.ratio <- U_sampnext + log_bark.dens(samp_next, samp_current, grad_samp_next, delta) - U_sampcurr - 
      log_bark.dens(samp_current, samp_next, grad_samp_curr, delta)
    if(log(runif(1)) <= mh.ratio)
    {
      samp.bark[i,] <- samp_next
      U_sampcurr <- U_sampnext
      accept <- accept + 1
    }
    else
    {
      samp.bark[i,] <- samp_current
    }
    samp_current <- samp.bark[i,]
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(acc_rate <- accept/iter)
  object <- list(samp.bark, acc_rate)
  return(object)
}

##  MYHMC samples

myhmc <- function(y, alpha, lambda, sigma2, iter, eps_hmc, L, start)
{
  nvar <- length(y)
  samp.hmc <- matrix(0, nrow = iter, ncol = nvar)
  wts_is_est <- numeric(length = iter)
  
  # starting value computations
  samp <- start
  samp.hmc[1,] <- samp
  proxval_curr <- prox_func(samp, lambda, y, sigma2, alpha)
  
  # weights calculation
  psi_val <- - log_pi(samp,y,sigma2,alpha)
  psi_lambda_val <- - log_pilambda(proxval_curr,samp,lambda,y,sigma2,alpha)
  wts_is_est[1] <- psi_lambda_val - psi_val
  
  # For HMC
  mom_mat <- matrix(rnorm(iter*nvar), nrow = iter, ncol = nvar)
  accept <- 0
  
  for (i in 2:iter) 
  {
    p_prop <- mom_mat[i,]
    U_samp <- -grad_logpiLam(samp, lambda,y,sigma2,alpha)
    p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
    q_current <- samp
    for (j in 1:L)
    {
      samp <- samp + eps_hmc*p_current   # full step for position
      U_samp <- -grad_logpiLam(samp, lambda,y,sigma2,alpha)
      if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
    }
    p_current <- p_current - eps_hmc*U_samp/2
    p_current <- - p_current  # negation to make proposal symmetric
    
    #  proximal values
    proxval_prop <- prox_func(samp,lambda,y,sigma2,alpha)
    U_curr <- - log_pilambda(proxval_curr,q_current,lambda,y,sigma2,alpha)
    U_prop <- - log_pilambda(proxval_prop,samp,lambda,y,sigma2,alpha)
    K_curr <-  sum((p_prop^2)/2)
    K_prop <-  sum((p_current^2)/2)
    
    log_acc_prob = U_curr - U_prop + K_curr - K_prop
    
    if(log(runif(1)) <= log_acc_prob )
    {
      samp.hmc[i,] <- samp
      proxval_curr <- proxval_prop
      
      # weights
      psi_val <- - log_pi(samp,y,sigma2,alpha)
      wts_is_est[i] <- U_prop - psi_val
      accept <- accept + 1
    }
    else
    {
      samp.hmc[i,] <- q_current
      psi_val <- - log_pi(q_current,y,sigma2,alpha)
      wts_is_est[i] <- U_curr - psi_val
      samp <- q_current
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

pxhmc <- function(y, alpha, lambda, sigma2, iter, eps_hmc, L, start)
{
  nvar <- length(y)
  samp.hmc <- matrix(0, nrow = iter, ncol = nvar)
  
  # starting value computations
  samp <- start
  samp.hmc[1,] <- samp
  
  # For HMC
  mom_mat <- matrix(rnorm(iter*length(y)), nrow = iter, ncol = length(y))
  accept <- 0
  
  for (i in 2:iter) 
  {
    p_prop <- mom_mat[i,]
    U_samp <- -grad_logpiLam(samp, lambda,y,sigma2,alpha)
    p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
    q_current <- samp
    for (j in 1:L)
    {
      samp <- samp + eps_hmc*p_current   # full step for position
      U_samp <- -grad_logpiLam(samp, lambda,y,sigma2,alpha)
      if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
    }
    p_current <- p_current - eps_hmc*U_samp/2
    p_current <- - p_current  # negation to make proposal symmetric
    
    U_curr <- sum((y - q_current)^2)/(2*sigma2) + alpha*nucl_norm(q_current)
    U_prop <- sum((y - samp)^2)/(2*sigma2) + alpha*nucl_norm(samp)
    K_curr <-  sum((p_prop^2)/2)
    K_prop <-  sum((p_current^2)/2)
    
    log_acc_prob = U_curr - U_prop + K_curr - K_prop
    
    if(log(runif(1)) <= log_acc_prob )
    {
      samp.hmc[i,] <- samp
      accept <- accept + 1
    }
    else
    {
      samp.hmc[i,] <- q_current
      samp <- q_current
    }
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))}
  } 
  print(acc_rate <- accept/iter)
  object <- list(samp.hmc, acc_rate)
  return(object)
}