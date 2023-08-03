
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
  
px.mala <- function(in_val, iter, lambda, delta)
{
  samp.pxm <- numeric(length = iter)
  samp.pxm[1] <- in_val
  
  for (i in 2:iter) 
  {
    propval <- rnorm(1, in_val + (delta / 2)*log_gradpi(in_val, lambda), sqrt(delta))   # proposal step
    targ_val.pr <- - prox_arg(propval, prox_func(propval, lambda), lambda)
    targ_val.in <- - prox_arg(in_val, prox_func(in_val, lambda), lambda)
    mh.ratio <- targ_val.pr + dnorm(in_val, propval + (delta / 2)*log_gradpi(propval, lambda), sqrt(delta), log = TRUE) -
      (targ_val.in + dnorm(propval, in_val + (delta / 2)*log_gradpi(in_val, lambda), sqrt(delta), log = TRUE))
    if(log(runif(1)) <= mh.ratio)
    {
      samp.pxm[i] <- propval
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
  
  for (i in 2:iter)
  {
    propval <- bark.prop(in_val, delta, lambda)
    targ_val.pr <- - prox_arg(propval, prox_func(propval, lambda), lambda)
    targ_val.in <- - prox_arg(in_val, prox_func(in_val, lambda), lambda)
    mh.ratio <- targ_val.pr + log(bark.dens(propval, in_val, delta, lambda)) - targ_val.in -
      log(bark.dens(in_val, propval, delta, lambda))
    if(log(runif(1)) <= mh.ratio)
    {
      samp.bark[i] <- propval
    }
    else
    {
      samp.bark[i] <- in_val
    }
    in_val <- samp.bark[i]
  }
  return(samp.bark)
}


iter <- 1e5
in_val <- 1
delta <- 1
lambda.vec <- c(1, 0.5, 0.005, 0.00001)

#  MYMALA samples 

mala.l1 <- px.mala(in_val, iter, lambda.vec[1], delta)
mala.l2 <- px.mala(in_val, iter, lambda.vec[2], delta)
mala.l3 <- px.mala(in_val, iter, lambda.vec[3], delta)
mala.l4 <- px.mala(in_val, iter, lambda.vec[4], delta)
mala <- list(mala.l1, mala.l2, mala.l3, mala.l4)

# MYBarker samples

barker.l1 <- px.barker(in_val, iter, lambda.vec[1], delta)
barker.l2 <- px.barker(in_val, iter, lambda.vec[2], delta)
barker.l3 <- px.barker(in_val, iter, lambda.vec[3], delta)
barker.l4 <- px.barker(in_val, iter, lambda.vec[4], delta)
barker <- list(barker.l1, barker.l2, barker.l3, barker.l4)

# Density plots

par(mfrow = c(2,2))

for (k in 1:4)
{
  plot(density(sample <- rlaplace(1e5, 0, 1)), xlab = "values",
       ylab = "density", main = paste("lambda = ", lambda.vec[k]), col = "black") 
  lines(density(as.numeric(unlist(mala[k]))) ,col = "blue")
  lines(density(as.numeric(unlist(barker[k]))) ,col = "red")
  legend("topright", c("Truth", "MYMALA", "MYBarker"),
         col = c("black", "blue", "red"), cex = 0.5,fill = c("black", "blue", "red"))
}

# acf plots

par(mfrow = c(2,2))

for (k in 1:4)
{
  acf.mala <- acf(as.numeric(unlist(mala[k])), plot = FALSE)$acf
  acf.barker <- acf(as.numeric(unlist(barker[k])), plot = FALSE)$acf
  plot(1:length(acf.mala), acf.mala, col = "blue", type = 'l',
       xlab = "Lag", ylab = "Autocorrelation", main = paste("ACF plot for lambda = ", lambda.vec[k]))
  lines(1:length(acf.mala), acf.barker, col = "red", type = 'l')
  legend("topright", c("Truth", "MYMALA", "MYBarker"),
         col = c("black", "blue", "red"), cex = 0.5,fill = c("black", "blue", "red"))
}

# Trace plots

par(mfrow = c(2,2))

for (k in 1:4)
{
  plot.ts(as.numeric(unlist(mala[k])), col = "blue", ylab = "sample values",
          main = paste("Trace plot for lambda = ", lambda.vec[k]))
  lines(as.numeric(unlist(barker[k])), col = "red")
  legend("topright", c("Truth", "MYMALA", "MYBarker"),
         col = c("black", "blue", "red"), cex = 0.5,fill = c("black", "blue", "red"))
}
