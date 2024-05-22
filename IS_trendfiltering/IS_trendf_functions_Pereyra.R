
## Code for asymptotic variance comparison of weighted importance sampling estimator 
## for trend filtering example
library(mcmcse)
library(glmgen)
library(Matrix)
library(expm)
library(foreach)
library(doParallel)

set.seed(12345)
alpha_hat <- 5   # obtained from the first dataset
sigma2_hat <- 9  # obtained from the first dataset
x <- seq(1,100,len=100)
 f <- Vectorize(function(x){if(x<=35){x} else if(x<=70){70-x} else{0.5*x-35}})
fx_linear <- f(x)
y <- fx_linear + rnorm(length(x), sd = 3)
tol <- 1e-6
tol_pcm <- 1e-6
max_iter <- 200L
max_iter_pcm <- 200L

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

bark.prop <- function(beta,lambda,y,sigma2,alpha,k,grid,delta,covmat)
{
  aux_var <- rnorm(length(beta), 0, 1)
  y <- (delta*covmat)%*%aux_var
  denom_prod <- y*log_gradpi(beta,lambda,y,sigma2,alpha,k,grid)
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


#### MYMALA sampling for covariance matrix estimation

mymala_cov_fn <- function(y, alpha, sigma2, k, grid, iter, delta) # pre-conditioned mala
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

mymala <- function(y, alpha, sigma2, k, grid, iter, delta, covmat)
{
  samp.mym <- matrix(0, nrow = iter, ncol = length(y))
  lambda <- lamb_coeff
  wts_is_est <- numeric(length = iter)
  beta_current <- markov_chain[1e5,]
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

px.mala <- function(y, alpha, sigma2, k, grid, iter, delta, covmat)
{
  samp.pxm <- matrix(0, nrow = iter, ncol = length(y))
  lambda <- lamb_coeff*sigma2
  beta_current <- markov_chain[1e5,]
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
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(accept/iter)
  return(samp.pxm)
}

## IS from Barker

mybarker <- function(y, alpha, sigma2, k, grid, iter, delta, covmat)
{
  samp.bark <- matrix(0, nrow = iter, ncol = length(y))
  lambda <- lamb_coeff
  wts_is_est <- numeric(length = iter)
  beta_current <- markov_chain[1e5,]
  samp.bark[1,] <- beta_current
  g_val <- alpha*sum(abs(D_mat%*%beta_current)) + sum((y - beta_current)^2)/(2*sigma2)
  prox_val.curr <- prox_func(beta_current, lambda, alpha, sigma2, k, grid)
  g_lambda_val <- prox_arg(prox_val.curr, beta_current, lambda=lambda, y, sigma2, alpha)
  wts_is_est[1] <- g_lambda_val - g_val
  U <- sqrtm(covmat)
  mat.inv <- solve(covmat)
  accept <- 0
  for (i in 2:iter) 
  {
    beta_next <- bark.prop(beta_current,lambda,y,sigma2,alpha,k,grid,delta,covmat) # proposal step
    prox_val.next <- prox_func(beta_next, lambda, alpha, sigma2, k, grid)
    targ_val.next <- log_target(prox_val.next,beta_next,lambda,y,sigma2,alpha)
    targ_val.curr <- log_target(prox_val.curr,beta_current,lambda,y,sigma2,alpha)
    grad_beta_curr <- log_gradpi(beta_current, lambda, y, sigma2, alpha, k, grid)
    grad_beta_next <- log_gradpi(beta_next, lambda, y, sigma2, alpha, k, grid)
    
    mh.ratio <- targ_val.next - targ_val.curr + 
      log_q_ratio_barker(beta_current,beta_next,grad_beta_curr,grad_beta_next)
    
    if(log(runif(1)) <= mh.ratio)
    {
      samp.bark[i,] <- beta_next
      prox_val.curr <- prox_val.next
      g_val <- alpha*sum(abs(D_mat%*%beta_next)) + sum((y - beta_next)^2)/(2*sigma2)
      g_lambda_val <- prox_arg(prox_val.next, beta_next, lambda=lambda, y, sigma2, alpha)
      wts_is_est[i] <- g_lambda_val - g_val
      accept <- accept + 1
    }
    else
    {
      samp.bark[i,] <- beta_current
      g_val <- alpha*sum(abs(D_mat%*%beta_current)) + sum((y - beta_current)^2)/(2*sigma2)
      g_lambda_val <- prox_arg(prox_val.curr, beta_current, lambda=lambda, y, sigma2, alpha)
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

px.barker <- function(y, alpha, sigma2, k, grid, iter, delta, covmat)
{
  samp.bark <- matrix(0, nrow = iter, ncol = length(y))
  lambda <- lamb_coeff
  beta_current <- markov_chain[1e5,]
  samp.bark[1,] <- beta_current
  mat.inv <- solve(covmat)
  accept <- 0
  for (i in 2:iter)
  {
    beta_next <- bark.prop(beta_current,lambda,y,sigma2,alpha,k,grid,delta,covmat)
    U_betanext <- - (sum((y - beta_next)^2)/(2*sigma2) + alpha*(sum(abs(D_mat%*%beta_next))))
    U_betacurr <- - (sum((y - beta_current)^2)/(2*sigma2) + alpha*(sum(abs(D_mat%*%beta_current))))
    grad_beta_curr <- log_gradpi(beta_current, lambda, y, sigma2, alpha, k, grid)
    grad_beta_next <- log_gradpi(beta_next, lambda, y, sigma2, alpha, k, grid)
    
    mh.ratio <- U_betanext - U_betacurr + 
      log_q_ratio_barker(beta_current,beta_next,grad_beta_curr,grad_beta_next)
    if(log(runif(1)) <= mh.ratio)
    {
      samp.bark[i,] <- beta_next
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

myhmc <- function(y, alpha, sigma2, k, grid, iter, eps_hmc, L)
{
  samp.hmc <- matrix(0, nrow = iter, ncol = length(y))
  lambda <- lamb_coeff
  wts_is_est <- numeric(length = iter)
  beta <- markov_chain[1e5,]
  samp.hmc[1,] <- beta
  g_val <- alpha*sum(abs(D_mat%*%beta)) + sum((y - beta)^2)/(2*sigma2)
  proxval_curr <- prox_func(beta, lambda, alpha, sigma2, k, grid)
  g_lambda_val <- prox_arg(proxval_curr, beta, lambda=lambda, y, sigma2, alpha)
  wts_is_est[1] <- g_lambda_val - g_val
  mom_mat <- matrix(rnorm(iter*length(y)), nrow = iter, ncol = length(y))
  accept <- 0
  
  for (i in 2:iter) 
  {
    p_prop <- mom_mat[i,]
    U_beta <- -log_gradpi(beta, lambda,y,sigma2,alpha,k,grid)
    p_current <- p_prop - eps_hmc*U_beta /2  # half step for momentum
    q_current <- beta
    for (j in 1:L)
    {
      beta <- beta + eps_hmc*p_current   # full step for position
      U_beta <- -log_gradpi(beta, lambda,y,sigma2,alpha,k,grid)
      if(j!=L) p_current <- p_current - eps_hmc*U_beta  # full step for momentum
    }
    p_current <- p_current - eps_hmc*U_beta/2
    p_current <- - p_current  # negation to make proposal symmetric
    
    #  proximal values
    proxval_prop <- prox_func(beta = beta,lambda,alpha,sigma2,k,grid)
    
    U_curr <- - log_target(proxval_curr,q_current,lambda,y,sigma2,alpha)
    U_prop <- - log_target(proxval_prop,beta,lambda,y,sigma2,alpha)
    K_curr <-  sum((p_prop^2)/2)
    K_prop <-  sum((p_current^2)/2)
    
    log_acc_prob = U_curr - U_prop + K_curr - K_prop
    
    if(log(runif(1)) <= log_acc_prob )
    {
      samp.hmc[i,] <- beta
      proxval_curr <- proxval_prop
      g_val <- alpha*sum(abs(D_mat%*%beta)) + sum((y - beta)^2)/(2*sigma2)
      g_lambda_val <- prox_arg(proxval_prop, beta, lambda=lambda, y, sigma2, alpha)
      wts_is_est[i] <- g_lambda_val - g_val
      accept <- accept + 1
    }
    else
    {
      samp.hmc[i,] <- q_current
      g_val <- alpha*sum(abs(D_mat%*%q_current)) + sum((y - q_current)^2)/(2*sigma2)
      g_lambda_val <- prox_arg(proxval_curr, q_current, lambda=lambda, y, sigma2, alpha)
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

pxhmc <- function(y, alpha, sigma2, k, grid, iter, eps_hmc, L)
{
  samp.hmc <- matrix(0, nrow = iter, ncol = length(y))
  lambda <- lamb_coeff
  beta <- markov_chain[1e5,]
  samp.hmc[1,] <- beta
  mom_mat <- matrix(rnorm(iter*length(y)), nrow = iter, ncol = length(y))
  accept <- 0
  
  for (i in 2:iter) 
  {
    p_prop <- mom_mat[i,]
    U_beta <- -log_gradpi(beta, lambda,y,sigma2,alpha,k,grid)
    p_current <- p_prop - eps_hmc*U_beta /2  # half step for momentum
    q_current <- beta
    for (j in 1:L)
    {
      beta <- beta + eps_hmc*p_current   # full step for position
      U_beta <- -log_gradpi(beta, lambda,y,sigma2,alpha,k,grid)
      if(j!=L) p_current <- p_current - eps_hmc*U_beta  # full step for momentum
    }
    p_current <- p_current - eps_hmc*U_beta/2
    p_current <- - p_current  # negation to make proposal symmetric
    
    U_curr <- sum((y - q_current)^2)/(2*sigma2) + alpha*(sum(abs(D_mat%*%q_current)))
    U_prop <- sum((y - beta)^2)/(2*sigma2) + alpha*(sum(abs(D_mat%*%beta)))
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
