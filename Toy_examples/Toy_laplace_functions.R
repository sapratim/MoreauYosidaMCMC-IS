
## Code for Px-MALA and Px-Barker algorithms with standard Laplace target

library(extraDistr)

# function calculates the inside of the proximal function

prox_arg <- function(x, y, mu)     # x is the current value
{
  z_x <- abs(y) + (x - y)^2 / (2*mu)    
  return(z_x)
}

# function calculates the value of the proximal function

prox_func <- function(val, mu)
{
  vec <- c(0, val - mu, val + mu)
  fun.val <- prox_arg(val, vec, mu)
  index <- which.min(fun.val)
  prox <- vec[index]
  return(prox)
}

grad_logpi <- function(val, mu)
{
  gradval <- - (val - prox_func(val, mu)) / mu 
  return(gradval)
}

# bark.prop <- function(val, delta, lambda)
# {
#   y <- rnorm(1, 0, sqrt(delta))
#   prob <- 1 / (1 + exp( -y*grad_logpi(val, lambda)))
#   ifelse(runif(1) <= prob, prop <- val + y, prop <- val - y)
#   return(prop)
# }
# 
# bark.dens <- function(in_val, propval, delta, lambda)
# {
#   numer <- 2 * dnorm(propval, in_val, sqrt(delta))
#   denom <- 1 + exp( - (propval - in_val) * grad_logpi(in_val, lambda))
#   value <- numer / denom
#   return(value)
# }
  
mymala <- function(curr_val, iter, lambda, delta)
{
  samp.mym <- numeric(length = iter)
  wts_is <- numeric(length = iter)
  
  samp.mym[1] <- curr_val
  proxval_curr <- prox_func(curr_val, lambda)
  targ_val.curr <- - prox_arg(curr_val, proxval_curr, lambda)
  wts_is[1] <-  prox_arg(curr_val, proxval_curr, lambda) - abs(curr_val) 
  
  accept <- 0
  for (i in 2:iter) 
  {
    prop.mean <- curr_val + (delta / 2)*grad_logpi(curr_val, lambda)
    propval <- prop.mean + sqrt(delta)*rnorm(1)   # proposal step
    
    proxval_next <- prox_func(propval, lambda)
    targ_val.next <- - prox_arg(propval, proxval_next, lambda)
    
    mh.ratio <- targ_val.next + dnorm(curr_val, propval + (delta / 2)*grad_logpi(propval, lambda), 
                                      sqrt(delta), log = TRUE) - targ_val.curr - 
      dnorm(propval, prop.mean, sqrt(delta), log = TRUE)
    if(log(runif(1)) <= mh.ratio)
    {
      samp.mym[i] <- propval
      accept <- accept + 1
      proxval_curr <- proxval_next
      targ_val.curr <- targ_val.next
      
      # weights
      wts_is[i] <- - targ_val.curr - abs(propval)
    }
    else
    {
      samp.mym[i] <- curr_val
      wts_is[i] <-  - targ_val.curr - abs(propval) 
    }
    curr_val <- samp.mym[i]
  }
  acceptance <- accept / iter
  print(paste("Acceptance rate of MY is = ", acceptance))
  output <- list(samp.mym, wts_is)
  return(output)
}

pxmala <- function(curr_val, iter, lambda, delta)
{
  samp.pxm <- numeric(length = iter)
  
  samp.pxm[1] <- curr_val
  targ_val.curr <- - abs(curr_val)
  
  accept <- 0
  for (i in 2:iter) 
  {
    prop.mean <- curr_val + (delta / 2)*grad_logpi(curr_val, lambda)
    propval <- prop.mean + sqrt(delta)*rnorm(1)   # proposal step
    
    targ_val.next <- - abs(propval)
    
    mh.ratio <- targ_val.next + dnorm(curr_val, propval + (delta / 2)*grad_logpi(propval, lambda), 
                                      sqrt(delta), log = TRUE) - targ_val.curr - 
      dnorm(propval, prop.mean, sqrt(delta), log = TRUE)
    if(log(runif(1)) <= mh.ratio)
    {
      samp.pxm[i] <- propval
      accept <- accept + 1
      targ_val.curr <- targ_val.next
    }
    else
    {
      samp.pxm[i] <- curr_val
    }
    curr_val <- samp.pxm[i]
  }
  acceptance <- accept / iter
  print(paste("Acceptance rate of Px is = ", acceptance))
  return(samp.pxm)
}


opt_lambda_func <- function(log_weights)
{
  true_wts <- exp(log_weights)
  ratio <- (mean(true_wts)^2)/mean(true_wts^2)
  return(ratio)
}

# to evaluate asymptotic covariance matrix
asymp_covmat_fn <- function(chain, weights)
{
  wts_mean <- mean(weights)
  num <- chain*weights
  is_est <- sum(num) / sum(weights)
  input_mat <- cbind(num, weights)  # input samples for mcse
  Sigma_mat <- mcse.multi(input_mat)$cov  # estimated covariance matrix of the tuple
  kappa_eta_mat <- cbind(1/wts_mean, -is_est/wts_mean) # derivative of kappa matrix
  asymp_covmat <- (kappa_eta_mat %*% Sigma_mat) %*% t(kappa_eta_mat) # IS asymptotic variance
  return(asymp_covmat)
}
