
## Code for weighted importance sampling estimator for standard Laplace target

set.seed(123)
library(extraDistr)

# function calculates the inside of the proximal function

target_val <- function(x)
{
  u <- dlaplace(x, 0, 1, log = TRUE)
  return(u)
}

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
  wts_is_est[1] <- exp(-abs(in_val)) / exp(-prox_arg(in_val, prox_func(in_val, lambda), lambda))
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
      wts_is_est[i] <- (exp(-abs(propval)) / exp(targ_val.pr))
      accept <- accept + 1
    }
    else
    {
      samp.mym[i] <- in_val
      wts_is_est[i] <- (exp(-abs(in_val)) / exp(targ_val.in))
    }
    in_val <- samp.mym[i]
  }
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
  return(samp.bark)
}


iter <- 1e4
in_val <- 2
delta <- 1
lambda.vec <- c(0.1, 1, 100, 500)

#  MYMALA samples and weights

mala.l1 <- mymala(in_val, iter, lambda.vec[1], delta)
mala.l2 <- mymala(in_val, iter, lambda.vec[2], delta)
mala.l3 <- mymala(in_val, iter, lambda.vec[3], delta)
mala.l4 <- mymala(in_val, iter, lambda.vec[4], delta)

# PXMALA samples

pxmala.l1 <- px.mala(in_val, iter, lambda.vec[1], delta)
pxmala.l2 <- px.mala(in_val, iter, lambda.vec[1], delta)
pxmala.l3 <- px.mala(in_val, iter, lambda.vec[1], delta)
pxmala.l4 <- px.mala(in_val, iter, lambda.vec[1], delta)

# PxBarker samples

barker.l1 <- px.barker(in_val, iter, lambda.vec[1], delta)
barker.l2 <- px.barker(in_val, iter, lambda.vec[1], delta)
barker.l3 <- px.barker(in_val, iter, lambda.vec[1], delta)
barker.l4 <- px.barker(in_val, iter, lambda.vec[1], delta)


mala_samp <- list(mala.l1[1], mala.l2[1], mala.l3[1], mala.l4[1])
mala_wts <- list(mala.l1[2], mala.l2[2], mala.l3[2], mala.l4[2])
pxm_samp <- list(pxmala.l1, pxmala.l2, pxmala.l3, pxmala.l4)
barker_samp <- list(barker.l1, barker.l2, barker.l3, barker.l4)

# sample size vector and mean matrix initialisation

samp <- c(10, 1e2, 1e3, 1e4)
mean_mat.is <- matrix(0, nrow = length(samp), ncol = length(lambda.vec))
mean_mat.pxm <- matrix(0, nrow = length(samp), ncol = length(lambda.vec))
mean_mat.bark <- matrix(0, nrow = length(samp), ncol = length(lambda.vec))


for (i in 1:length(lambda.vec)) 
  {
   sam_is <- as.numeric(unlist(mala_samp[i]))
   weights <- as.numeric(unlist(mala_wts[i]))
   sam_pxm <- as.numeric(unlist(pxm_samp[i]))
   sam_bark <- as.numeric(unlist(barker_samp[i]))
   
   for (j in 1:length(samp))
    {
       num <- sum(abs(sam_is[1:samp[j]])*weights[1:samp[j]])   
       mean_mat.is[j, i] <- num / (sum(weights[1:samp[j]])) 
       mean_mat.pxm[j, i] <- mean(sam_pxm[1:samp[j]])
       mean_mat.bark[j, i] <- mean(sam_bark[1:samp[j]])
   }  
}      

mean_mat.is      
mean_mat.pxm
mean_mat.bark
