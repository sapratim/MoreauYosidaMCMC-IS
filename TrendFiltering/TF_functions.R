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

# function for barker proposal
bark.prop <- function(beta,lambda,y,sigma2,alpha,k,grid,delta)
{
  aux_var <- rnorm(length(beta), 0, 1)
  z <- sqrt(delta)*aux_var
  denom_prod <- z*grad_logpiLam(beta,lambda,y,sigma2,alpha,k,grid)
  prob <- 1 / (1 + exp(- denom_prod))
  unifs <- runif(length(beta))
  prop <- (beta + z)*(unifs <= prob) + (beta - z)*(unifs > prob)
  return(prop)
}

# barker log density
log_bark.dens <- function(curr_point, prop_point, grad_curr_point, delta)
{
  rw_dens <- prod(dnorm(prop_point, curr_point, sqrt(delta)))
  exp_term <- - (grad_curr_point*(prop_point-curr_point))
  denom <- prod(1 + exp(exp_term))
  dens_val <- rw_dens/denom
  return(log(dens_val))
}

# to evaluate asymptotic covariance matrix
asymp_covmat_fn <- function(chain, weights)
{
  wts_mean <- mean(weights)
  num <- chain*weights
  sum_mat <- apply(num, 2, sum)
  is_est <- sum_mat / sum(weights)
  input_mat <- cbind(num, weights)  # input samples for mcse
  Sigma_mat <- mcse.multi(input_mat)$cov  # estimated covariance matrix of the tuple
  kappa_eta_mat <- cbind(diag(1/wts_mean, ncol(chain)), -is_est/wts_mean) # derivative of kappa matrix
  asymp_covmat <- (kappa_eta_mat %*% Sigma_mat) %*% t(kappa_eta_mat) # IS asymptotic variance
  return(asymp_covmat)
}

# functions for quantile estimation
quant <- function(j, mat)     
{
  mat_ordered <- mat[order(mat[,j], decreasing = FALSE), ]
  order_comp <- mat_ordered[,j]
  weights_order <- mat_ordered[, length(y)+1]
  bound_mat <- cbind(order_comp, weights_order)
  return(bound_mat)
}

quantile_func <- function(chain, weights, signif_level)
{
  augm_mat <- cbind(chain,weights)
  upper_quant <- numeric(length = ncol(chain))
  lower_quant <- numeric(length = ncol(chain))
  post_med <- numeric(length = ncol(chain))

  for (i in 1:length(y)) 
  {
   initial_mat <- quant(i, augm_mat)
   mat_sum <- apply(initial_mat, 2, sum)
   wts_prop <- initial_mat[,2]/mat_sum[2]
   final_mat <- cbind(initial_mat[1,], cumsum(wts_prop))
   lower_index <- min(which(final_mat[,2] >= signif_level))
   upper_index <- min(which(final_mat[,2] >= (1 - signif_level)))
   med_index <- min(which(final_mat[,2] >= 0.5))
   upper_quant[i] <- initial_mat[upper_index,1]
   lower_quant[i] <- initial_mat[lower_index,1]
   post_med[i] <- initial_mat[med_index,1]
  }
  quantiles <- list(upper_quant, lower_quant, post_med)
  return(quantiles)
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
  
  # for MALA
  accept <- 0
  for (i in 2:iter) 
  {
    # proposal step
    prop.mean <- beta_current + (delta / 2)*grad_logpiLam(beta_current,lambda,y,sigma2,alpha,k,grid)
    beta_next <-  prop.mean + (sqrt(delta))*rnorm(nvar, 0, 1)  
    
    # calculating prox
    prox_val.next <- prox_func(beta_next, lambda, alpha, sigma2, k, grid)
    targ_val.next <- log_pilambda(prox_val.next,beta_next,lambda,y,sigma2,alpha)
    
    q.next_to_curr <- sum(dnorm(beta_current, beta_next + 
                                   (delta / 2)*grad_logpiLam(beta_next,lambda,y,sigma2,alpha,k,grid),
                                  sqrt(delta), log = TRUE))
    q.curr_to_next <- sum(dnorm(beta_next, prop.mean, sqrt(delta), log = TRUE)) 
    
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
  
  for (i in 2:iter)
  {
    # proposal step
    prop.mean <- beta_current +  (delta / 
                                    2)*grad_logpiLam(beta_current,lambda,y,sigma2,alpha,k,grid)
    beta_next <-  prop.mean + sqrt(delta)*rnorm(nvar, 0, 1) 
    
    U_betanext <- log_pi(beta_next, y, sigma2, alpha)
    q.next_to_curr <- sum(dnorm(beta_current, beta_next +
                                   (delta /  2)*grad_logpiLam
                                 (beta_next,lambda,y,sigma2,alpha,k,grid), 
                                 sqrt(delta), log = TRUE))
    q.curr_to_next <- sum(dnorm(beta_next, prop.mean, sqrt(delta), log = TRUE))
    
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
  
  # For barker  
  for (i in 2:iter) 
  {
    # proposal step
    beta_next <- bark.prop(beta_current,lambda,y,sigma2,alpha,k,grid,delta)
    
    prox_val.next <- prox_func(beta_next, lambda, alpha, sigma2, k, grid)
    targ_val.next <- log_pilambda(prox_val.next,beta_next,lambda,y,sigma2,alpha)
    
    grad_beta_curr <- grad_logpiLam(beta_current, lambda, y, sigma2, alpha, k, grid)
    grad_beta_next <- grad_logpiLam(beta_next, lambda, y, sigma2, alpha, k, grid)
    
    mh.ratio <- targ_val.next + log_bark.dens(beta_next, beta_current, grad_beta_next, delta) - targ_val.curr - 
               log_bark.dens(beta_current, beta_next, grad_beta_curr, delta)
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
  samp.bark <- matrix(0, nrow = iter, ncol = nvar)
  lambda <- lamb_coeff
  
  # starting value computations
  beta_current <- start
  samp.bark[1,] <- beta_current
  U_betacurr <- log_pi(beta_current, y, sigma2, alpha)
  
  # For barker
  accept <- 0
  for (i in 2:iter)
  {
    # proposal step
    beta_next <- bark.prop(beta_current,lambda,y,sigma2,alpha,k,grid,delta)
    
    U_betanext <- log_pi(beta_next, y, sigma2, alpha)
    grad_beta_curr <- grad_logpiLam(beta_current, lambda, y, sigma2, alpha, k, grid)
    grad_beta_next <- grad_logpiLam(beta_next, lambda, y, sigma2, alpha, k, grid)
    
    mh.ratio <- U_betanext + log_bark.dens(beta_next, beta_current, grad_beta_next, delta) - U_betacurr - 
              log_bark.dens(beta_current, beta_next, grad_beta_curr, delta)
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
  
  # For HMC
  mom_mat <- matrix(rnorm(iter*nvar), nrow = iter, ncol = nvar)
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
      wts_is_est[i] <- U_prop - g_val
      accept <- accept + 1
    }
    else
    {
      samp.hmc[i,] <- q_current
      g_val <- -log_pi(q_current, y, sigma2, alpha)
      wts_is_est[i] <- U_curr - g_val
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
  
  # For HMC
  mom_mat <- matrix(rnorm(iter*nvar), nrow = iter, ncol = nvar)
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














##### MYMALA samples function

mymala <- function(eta_start, mu_start, lambda, sigma, iter, delta, data)
{
  samp.mym <- matrix(0, nrow = iter, ncol = I+1)
  wts_is_est <- numeric(length = iter)
  
  # starting values
  samp_current <- c(eta_start, mu_start)
  samp.mym[1,] <- samp_current
  prox_val.curr <- proxfunc(eta_start, mu_start, lambda, eta_start, mu_start, sigma)
  targ_val.curr <- log_plam(eta_start, mu_start, data, lambda,
                           c(prox_val.curr[1:I], prox_val.curr[I+1]))
  
  # weights calculation
  psi_val <- - log_p(eta_start, mu_start, data)
  psi_lambda_val <- - targ_val.curr
  wts_is_est[1] <-  psi_lambda_val -  psi_val
  # psi_lambda_val <- log_p(prox_val.curr[1:I], prox_val.curr[I+1], data) -
  #   sum((prox_val.curr - c(eta_start, mu_start))^2)/(2*lambda)
  
   accept <- 0
  for (i in 2:iter)
  {
    # proposal step
     prop.mean <- samp_current +
           (delta / 2)*grad_logplam(samp_current[1:I], samp_current[I+1],
                                lambda, eta_start, mu_start, sigma)
    samp_next <- rnorm(length(samp_current), prop.mean,  sqrt(delta))
    
    # calculating prox values
    prox_val.next <- proxfunc(samp_next[1:I], samp_next[I+1],
                               lambda, eta_start, mu_start, sigma)
    
    targ_val.next <- log_plam(samp_next[1:I], samp_next[I+1], data, lambda,
                              c(prox_val.next[1:I], prox_val.next[I+1]))
    
    q.next_to_curr <- sum(dnorm(samp_current, samp_next +
                       (delta / 2)*grad_logplam(samp_next[1:I], samp_next[I+1],
                         lambda, eta_start, mu_start, sigma),sqrt(delta), log = TRUE))
   
     q.curr_to_next <- sum(dnorm(samp_next, prop.mean,  sqrt(delta), log = TRUE))
    
    mh.ratio <- targ_val.next + q.next_to_curr - (targ_val.curr + q.curr_to_next)  # mh  ratio
    
    if(log(runif(1)) <= mh.ratio)
    {
      samp.mym[i,] <- samp_next
      prox_val.curr <- prox_val.next
      targ_val.curr <- targ_val.next
      
      # weights
      psi_val <- - log_p(samp_next[1:I], samp_next[I+1], data)
      psi_lambda_val <- - targ_val.curr
      wts_is_est[i] <- psi_lambda_val - psi_val
      accept <- accept + 1
    }
    else
    {
      samp.mym[i,] <- samp_current
      wts_is_est[i] <- wts_is_est[i-1]
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

px.mala <- function(eta_start, mu_start, lambda, sigma, iter, delta, data)
{
  samp.pxm <- matrix(0, nrow = iter, ncol = I+1)
  
  #  starting values
  samp_current <- c(eta_start, mu_start)
  samp.pxm[1,] <- samp_current
  U_sampcurr <- log_p(samp_current[1:I], samp_current[I+1], data)
  accept <- 0
  
  for (i in 2:iter)
  {
    # proposal step
    prop.mean <- samp_current +
           (delta / 2)*grad_logplam(samp_current[1:I], samp_current[I+1],
                               lambda, eta_start, mu_start, sigma)
    samp_next <- rnorm(length(samp_current), prop.mean,  sqrt(delta))
    
    U_sampnext <- log_p(samp_next[1:I], samp_next[I+1], data)
    
    q.next_to_curr <- sum(dnorm(samp_current, samp_next +
                       (delta / 2)*grad_logplam(samp_next[1:I], samp_next[I+1],
                         lambda, eta_start, mu_start, sigma),sqrt(delta), log = TRUE))
    
    q.curr_to_next <- sum(dnorm(samp_next, prop.mean,  sqrt(delta), log = TRUE))
   
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

