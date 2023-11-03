
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
  print(accept/iter)
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
  print(accept/iter)
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
  print(accept/iter)
  return(samp.bark)
}

iter <- 1e3
in_val <- 2
lambda.vec <- c(0.1, 1, 100, 500)
delta_is <- c(2.8, 4, 255, 1275)
delta_pxm <- c(2.8, 3.2, 1, 1)
delta_bark <- c(5, 4.8, 1, 1)

# Sample variance 

reps <- 1e2
augm_mat_is <- matrix(0, nrow = reps, ncol = length(lambda.vec))
augm_mat_pxm <- matrix(0, nrow = reps, ncol = length(lambda.vec))
augm_mat_pxb <- matrix(0, nrow = reps, ncol = length(lambda.vec))

for (i in 1:length(lambda.vec)) 
  {
   for (j in 1:reps)
     {
     mala.is <- mymala(in_val, iter, lambda.vec[i], delta_is[i])
     px_mala <- px.mala(in_val, iter, lambda.vec[i], delta_pxm[i])
     px_bark <- px.barker(in_val, iter, lambda.vec[i], delta_bark[i])
     
     is_samp <- as.numeric(unlist(mala.is[1]))
     is_wts <- as.numeric(unlist(mala.is[2]))
     num <- sum(is_samp*is_wts)
     augm_mat_is[j, i] <- num / sum(is_wts) 
     augm_mat_pxm[j, i] <- mean(as.numeric(unlist(px_mala)))
     augm_mat_pxb[j, i] <- mean(as.numeric(unlist(px_bark)))
     }
  }

# matrix of second moments

sqmat_is <- augm_mat_is^2
sqmat_pxm <- augm_mat_pxm^2
sqmat_pxb <- augm_mat_pxb^2

# sum of matrix of second moments

secmom_mat_is <- apply(sqmat_is, 2, sum)
secmom_mat_pxm <- apply(sqmat_pxm, 2, sum)
secmom_mat_pxb <- apply(sqmat_pxb, 2, sum)

# mean square matrix of estimates

meansq_mat_is <- colMeans(augm_mat_is)^2
meansq_mat_pxm <- colMeans(augm_mat_pxm)^2
meansq_mat_pxb <- colMeans(augm_mat_pxb)^2

# sample variances

samp_var.is <- (secmom_mat_is - reps*meansq_mat_is) / (reps - 1)
samp_var.pxm <- (secmom_mat_pxm - reps*meansq_mat_pxm) / (reps - 1)
samp_var.pxb <- (secmom_mat_pxb - reps*meansq_mat_pxb) / (reps - 1)

var_mat <- rbind(samp_var.is, samp_var.pxm, samp_var.pxb)
colnames(var_mat) <- c("lambda = 0.1", "lambda = 1", "lambda = 100", "lambda = 500")

var_mat   # variance comparison

pdf("density_g_lambda_lap.pdf", height = 6, width = 12)
par(mfrow = c(2,2))

for (k in 1:4) {
  
  s <- mymala(in_val, iter, lambda.vec[k], delta_is[k])
  u <- as.numeric(unlist(s[1]))
  plot(density(u), main = bquote(lambda == .(lambda.vec[k])))
  
}
dev.off()
