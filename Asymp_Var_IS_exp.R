#set.seed(123)

library(extraDistr)
library(mcmcse)

target_val <- function(x)
{
  norm_const <- 2 * gamma(5/4)
  u <- - x^4 - log(norm_const)
  return(u)
}
# function calculates the inside of the proximal function

prox_arg <- function(x, y, mu)     # x is the current value
{
  z_x <-  y^4 + ((x-y)^2)/(2*mu)  
  return(z_x)
}

prox_func <- function(val, mu)
{
  t <- sqrt(3)*sqrt((mu^3)*(27*mu*(val^2) + 1)) + 9*(mu^2)*val
  numer <- (3^(1/3))*(t^(2/3)) - (3^(2/3))*mu
  denom <- 6*mu*(t^(1/3))
  prox <-  numer / denom   
  return(prox)
}

log_gradpi <- function(val, mu)
{
  gradval <- - (val - prox_func(val, mu)) / mu 
  return(gradval)
}

bark.prop <- function(val, delta, lambda)
{
  y <- rnorm(1, 0, sqrt(delta))
  prob <- 1 / (1 + exp( -y*log_gradpi(val, lambda)))
  ifelse(runif(1) <= prob, prop <- val + y, prop <- val - y)
  return(prop)
}

bark.dens <- function(in_val, propval, delta, lambda)
{
  numer <- 2 * dnorm(propval, in_val, sqrt(delta))
  denom <- 1 + exp( - (propval - in_val) * log_gradpi(in_val, lambda))
  value <- numer / denom
  return(value)
}

# MYMALA sampling and importance sampling weights

mymala <- function(in_val, iter, lambda, delta)
{
  samp.mym <- numeric(length = iter)
  wts_is_est <- numeric(length = iter)
  samp.mym[1] <- in_val
  wts_is_est[1] <- exp(-(in_val^4)) / exp(-prox_arg(in_val, prox_func(in_val, lambda), lambda))
  accept <- 0
  for (i in 2:iter) 
  {
    propval <- rnorm(1, in_val + (delta / 2)*log_gradpi(in_val, lambda), sqrt(delta))   # proposal step
    targ_val.pr <- - prox_arg(propval, prox_func(propval, lambda), lambda)
    targ_val.in <- - prox_arg(in_val, prox_func(in_val, lambda), lambda)
    mh.ratio <- targ_val.pr + dnorm(in_val, propval + (delta / 2)*log_gradpi(propval, lambda), sqrt(delta), log = TRUE) -
      (targ_val.in + dnorm(propval, in_val + (delta / 2)*log_gradpi(in_val, lambda), sqrt(delta), log = TRUE))
    if(log(runif(1)) <= mh.ratio)
    {
      samp.mym[i] <- propval
      wts_is_est[i] <- (exp(-(propval^4)) / exp(targ_val.pr))
      accept <- accept + 1
    }
    else
    {
      samp.mym[i] <- in_val
      wts_is_est[i] <- (exp(-(in_val^4)) / exp(targ_val.in))
    }
    in_val <- samp.mym[i]
  }
  #print(accept/iter)
  object <- list(samp.mym, wts_is_est)
  return(object)
}

px.mala <- function(in_val, iter, lambda, delta)
{
  samp.pxm <- numeric(length = iter)
  samp.pxm[1] <- in_val
  accept <- 0
  for (i in 2:iter) 
  {
    propval <- rnorm(1, in_val + (delta / 2)*log_gradpi(in_val, lambda), sqrt(delta))   # proposal step
    mh.ratio <- target_val(propval) + dnorm(in_val, propval + (delta / 2)*log_gradpi(propval, lambda), sqrt(delta), log = TRUE) - 
      target_val(in_val) - dnorm(propval, in_val + (delta / 2)*log_gradpi(in_val, lambda), sqrt(delta), log = TRUE)
    if(log(runif(1)) <= mh.ratio)
    {
      samp.pxm[i] <- propval
      accept <- accept + 1
    }
    else
    {
      samp.pxm[i] <- in_val
    }
    in_val <- samp.pxm[i]
  }
  #print(accept/iter)
  return(samp.pxm)
}

px.barker <- function(in_val, iter, lambda, delta)
{
  samp.bark <- numeric(length = iter)
  samp.bark[1] <- in_val
  accept <- 0
  for (i in 2:iter)
  {
    propval <- bark.prop(in_val, delta, lambda)
    mh.ratio <- target_val(propval) + log(bark.dens(propval, in_val, delta, lambda)) - target_val(in_val) -
      log(bark.dens(in_val, propval, delta, lambda))
    if(log(runif(1)) <= mh.ratio)
    {
      samp.bark[i] <- propval
      accept <- accept + 1
    }
    else
    {
      samp.bark[i] <- in_val
    }
    in_val <- samp.bark[i]
  }
  #print(accept/iter)
  return(samp.bark)
}

iter <- 1e4
in_val <- 2
lambda.vec <- c(0.1, 1, 100, 500)
delta_is <- c(1.2, 4.2, 260, 1300)
delta_pxm <- c(1, 0.7, 0.5, 0.45)
delta_bark <- c(1.2, 0.65, 0.5, 0.5)

asymp_covmat_is <- numeric(length = length(lambda.vec))
asymp_covmat_pxm <- numeric(length = length(lambda.vec))
asymp_covmat_pxb <- numeric(length = length(lambda.vec))

for (i in 1:length(lambda.vec))
{
  mala.is <- mymala(in_val, iter, lambda.vec[i], delta_is[i])
  px_mala <- px.mala(in_val, iter, lambda.vec[i], delta_pxm[i])
  px_bark <- px.barker(in_val, iter, lambda.vec[i], delta_bark[i])
  
  # Importance sampling estimator asymptotic variance
  
  is_samp <- as.numeric(unlist(mala.is[1]))
  is_wts <- as.numeric(unlist(mala.is[2]))
  wts_mean <- mean(is_wts)
  num <- is_samp*is_wts
  is_est <- sum(num) / sum(is_wts)
  input_mat <- cbind(num, is_wts)  # input samples for mcse
  Sigma_mat <- mcse.multi(input_mat)$cov  # estimated covarince matrix of the tuple
  kappa_eta_mat <- matrix(c(1 / wts_mean , is_est / wts_mean), 
                          nrow = 1, ncol = 2, byrow = TRUE)      # derivative of kappa matrix 
  asymp_covmat_is[i] <- (kappa_eta_mat %*% Sigma_mat) %*% t(kappa_eta_mat)
  
  # PxMALA asymptotic variance
  
  asymp_covmat_pxm[i] <- mcse.multi(px_mala)$cov
  
  # PxBarker asymptotic variance
  
  asymp_covmat_pxb[i] <- mcse.multi(px_bark)$cov 
}

# Asymptotic variance comparison

var_mat <- rbind(asymp_covmat_is, asymp_covmat_pxm, asymp_covmat_pxb)
colnames(var_mat) <- c("lambda = 0.1", "lambda = 1", "lambda = 100", "lambda = 500")

var_mat
