################################################################################
########################### Anisotropic Laplace ################################
################################################################################

log_pi <- function(beta, x)
{
  value <- sum(beta*abs(x))
  return(-value)
}

# log target of pi-lambda
log_pilambda <- function(eta,beta,x)
{
  dens_val <- sum(beta*abs(eta)) + sum((eta-x)^2)/(2*lambda)
  return(-dens_val)
}

#########  Soft threshold function

softthreshold <- function(beta, u, lambda) {       ####  u is a vector
  return(sign(u)* pmax(abs(u)-beta*lambda,0))
}

# function calculates the value of the proximal function
prox_func <- function(beta, x, lambda)
{
  prox_val <- softthreshold(beta, x, lambda)
  return(prox_val)
}

# gradient of log target (pi-lambda)
grad_logpiLam <- function(beta, x, lambda)  
{
  prox_val <- prox_func(beta, x, lambda)
  ans <-  (x-prox_val)/lambda
  return(-ans)
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

############### mymala #################

mymala <- function(beta, start, lambda, iter, delta)
{
  nvar <- length(start)
  samp.mym <- matrix(0, nrow = iter, ncol = nvar)
  wts_is_est <- numeric(length = iter)
  
  # starting value computations
  curr_val <- start
  prox_val.curr <- prox_func(beta, curr_val, lambda)
  targ_val.curr <- log_pilambda(prox_val.curr,beta,curr_val)
  samp.mym[1,] <- curr_val
  accept <- 0
  
  # weights calculation
  g_val <- -log_pi(beta, curr_val) 
  g_lambda_val <- -targ_val.curr
  wts_is_est[1] <- g_lambda_val - g_val
  
  # for MALA
  for (i in 2:iter) 
  {
    # proposal step
    prop.mean <- curr_val + (delta / 2)*grad_logpiLam(beta, curr_val, lambda)
    propval <-  prop.mean + (sqrt(delta))*rnorm(nvar, 0, 1)  
    
    # calculating prox
    prox_val.next <- prox_func(beta, propval, lambda)
    targ_val.next <- log_pilambda(prox_val.next,beta,propval)
    
    q.next_to_curr <- sum(dnorm(curr_val, propval + 
                                  (delta / 2)*grad_logpiLam(beta, propval, lambda),
                                sqrt(delta), log = TRUE))
    q.curr_to_next <- sum(dnorm(propval, prop.mean, sqrt(delta), log = TRUE)) 
    
    mh.ratio <- targ_val.next + q.next_to_curr - (targ_val.curr + q.curr_to_next)  # mh  ratio
   
    if(log(runif(1)) <= mh.ratio)
    {
      samp.mym[i,] <- propval
      targ_val.curr <- targ_val.next
      prox_val.curr <- prox_val.next
      
      # weights
      g_val <- -log_pi(beta, propval) 
      g_lambda_val <- - targ_val.curr
      wts_is_est[i] <- g_lambda_val - g_val
      accept <- accept + 1
    }
    else
    {
      samp.mym[i,] <- curr_val
      g_val <- -log_pi(beta, curr_val) 
      g_lambda_val <- - targ_val.curr
      wts_is_est[i] <- g_lambda_val - g_val
    }
    curr_val <- samp.mym[i,]
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(accept/iter)
  object <- list(samp.mym, wts_is_est)
  return(object)
}

############### pxmala #################

pxmala <- function(beta, start, lambda, iter, delta)
{
  nvar <- length(start)
  samp.pxm <- matrix(0, nrow = iter, ncol = nvar)
  
  # starting value computations
  curr_val <- start
  targ_val.curr <- log_pi(beta,curr_val)
  samp.pxm[1,] <- curr_val
  accept <- 0
  
  # for MALA
  for (i in 2:iter) 
  {
    # proposal step
    prop.mean <- curr_val + (delta / 2)*grad_logpiLam(beta, curr_val, lambda)
    propval <-  prop.mean + (sqrt(delta))*rnorm(nvar, 0, 1)  
    
    targ_val.next <- log_pi(beta,propval)
    
    q.next_to_curr <- sum(dnorm(curr_val, propval + 
                                  (delta / 2)*grad_logpiLam(beta, propval, lambda),
                                sqrt(delta), log = TRUE))
    q.curr_to_next <- sum(dnorm(propval, prop.mean, sqrt(delta), log = TRUE)) 
    
    mh.ratio <- targ_val.next + q.next_to_curr - (targ_val.curr + q.curr_to_next)  # mh  ratio
    
    if(log(runif(1)) <= mh.ratio)
    {
      samp.pxm[i,] <- propval
      targ_val.curr <- targ_val.next
      accept <- accept + 1
    }
    else
    {
      samp.pxm[i,] <- curr_val
    }
    curr_val <- samp.pxm[i,]
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(accept/iter)
  return(samp.pxm)
}

dimen_func <- function(d, lambda, iter, delta_is, delta_px)
{
  beta <- seq(1, d, by = 1)
  start <- rep(0.5, d)
  
  mym.output <- mymala(beta, start, lambda, iter, delta_is)
  pxm.output <- pxmala(beta, start, lambda, iter, delta_px)
  output <- list(mym.output, pxm.output)
  return(output)
}

############### myhmc #################

# myhmc <- function(y, alpha, sigma2, k, grid, iter, eps_hmc, L, start)
# {
#   nvar <- length(y)
#   samp.hmc <- matrix(0, nrow = iter, ncol = nvar)
#   lambda <- lamb_coeff
#   wts_is_est <- numeric(length = iter)
#   
#   # starting value computations
#   beta <- start
#   proxval_curr <- prox_func(beta, lambda, alpha, sigma2, k, grid)
#   samp.hmc[1,] <- beta
#   
#   # weights calculation
#   g_val <- - log_pi(beta, y, sigma2, alpha)
#   g_lambda_val <- - log_pilambda(proxval_curr, beta, lambda=lambda, y, sigma2, alpha)
#   wts_is_est[1] <- g_lambda_val - g_val
#   
#   # For HMC
#   mom_mat <- matrix(rnorm(iter*nvar), nrow = iter, ncol = nvar)
#   accept <- 0
#   for (i in 2:iter) 
#   {
#     p_prop <- mom_mat[i,]
#     U_beta <- -grad_logpiLam(beta, lambda,y,sigma2,alpha,k,grid)
#     p_current <- p_prop - eps_hmc*U_beta /2  # half step for momentum
#     q_current <- beta
#     for (j in 1:L)
#     {
#       beta <- beta + eps_hmc*p_current   # full step for position
#       U_beta <- -grad_logpiLam(beta, lambda,y,sigma2,alpha,k,grid)
#       if(j!=L) p_current <- p_current - eps_hmc*U_beta  # full step for momentum
#     }
#     p_current <- p_current - eps_hmc*U_beta/2
#     p_current <- - p_current  # negation to make proposal symmetric
#     
#     #  calculating prox
#     proxval_prop <- prox_func(beta = beta,lambda,alpha,sigma2,k,grid)
#     U_curr <- - log_pilambda(proxval_curr,q_current,lambda,y,sigma2,alpha)
#     U_prop <- - log_pilambda(proxval_prop,beta,lambda,y,sigma2,alpha)
#     K_curr <-  sum((p_prop^2)/2)
#     K_prop <-  sum((p_current^2)/2)
#     
#     log_acc_prob = U_curr - U_prop + K_curr - K_prop    # mh ratio
#     
#     if(log(runif(1)) <= log_acc_prob )
#     {
#       samp.hmc[i,] <- beta
#       proxval_curr <- proxval_prop
#       
#       # weights
#       g_val <- -log_pi(beta, y, sigma2, alpha)
#       wts_is_est[i] <- U_prop - g_val
#       accept <- accept + 1
#     }
#     else
#     {
#       samp.hmc[i,] <- q_current
#       g_val <- -log_pi(q_current, y, sigma2, alpha)
#       wts_is_est[i] <- U_curr - g_val
#       beta <- q_current
#     }
#     if(i %% (iter/10) == 0){
#       j <- accept/iter
#       print(cat(i, j))}
#   } 
#   print(acc_rate <- accept/iter)
#   object <- list(samp.hmc, wts_is_est, acc_rate)
#   return(object)
# }
# 
# ############### pxhmc #################
# 
# pxhmc <- function(y, alpha, sigma2, k, grid, iter, eps_hmc, L, start)
# {
#   nvar <- length(y)
#   samp.hmc <- matrix(0, nrow = iter, ncol = nvar)
#   lambda <- lamb_coeff
#   
#   # starting value computations
#   beta <- start
#   samp.hmc[1,] <- beta
#   
#   # For HMC
#   mom_mat <- matrix(rnorm(iter*nvar), nrow = iter, ncol = nvar)
#   accept <- 0
#   
#   for (i in 2:iter) 
#   {
#     p_prop <- mom_mat[i,]
#     U_beta <- -grad_logpiLam(beta, lambda,y,sigma2,alpha,k,grid)
#     p_current <- p_prop - eps_hmc*U_beta /2  # half step for momentum
#     q_current <- beta
#     for (j in 1:L)
#     {
#       beta <- beta + eps_hmc*p_current   # full step for position
#       U_beta <- -grad_logpiLam(beta, lambda,y,sigma2,alpha,k,grid)
#       if(j!=L) p_current <- p_current - eps_hmc*U_beta  # full step for momentum
#     }
#     p_current <- p_current - eps_hmc*U_beta/2
#     p_current <- - p_current  # negation to make proposal symmetric
#     
#     U_curr <- - log_pi(q_current, y, sigma2, alpha)
#     U_prop <- - log_pi(beta, y, sigma2, alpha)
#     K_curr <-  sum((p_prop^2)/2)
#     K_prop <-  sum((p_current^2)/2)
#     
#     log_acc_prob = U_curr - U_prop + K_curr - K_prop
#     
#     if(log(runif(1)) <= log_acc_prob )
#     {
#       samp.hmc[i,] <- beta
#       accept <- accept + 1
#     }
#     else
#     {
#       samp.hmc[i,] <- q_current
#       beta <- q_current
#     }
#     if(i %% (iter/10) == 0){
#       j <- accept/iter
#       print(cat(i, j))}
#   } 
#   print(acc_rate <- accept/iter)
#   object <- list(samp.hmc, acc_rate)
#   return(object)
# }
