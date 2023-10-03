#set.seed(123)

library(extraDistr)

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
  print(accept/iter)
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

iter <- 1e3
in_val <- 2
lambda.vec <- c(0.1, 1, 100, 500)
delta_is <- c(1.3, 4.2, 250, 1200)
delta_pxm <- c(1, 0.7, 0.5, 0.45)
delta_bark <- c(1.2, 0.65, 0.5, 0.5)

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

pdf("density_g_lambda.pdf", height = 6, width = 12)
par(mfrow = c(2,2))

for (k in 1:4) {
  
s <- mymala(in_val, 1e5, lambda.vec[k], delta_is[k])
u <- as.numeric(unlist(s[1]))
plot(density(u), main = bquote(lambda == .(lambda.vec[k])))

}
dev.off()