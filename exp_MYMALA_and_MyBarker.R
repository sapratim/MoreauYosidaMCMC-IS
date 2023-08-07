## Code for Px-MALA and Px-Barker algorithms with standard Laplace target

# function calculates the inside of the proximal function

prox_arg <- function(x, y, mu)     # x is the current value
{
  z_x <-  y^4 + ((x-y)^2)/2*mu  
  return(z_x)
}

# function calculates the value of the proximal function

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

mymala <- function(in_val, iter, lambda, delta)
{
  samp.pxm <- numeric(length = iter)
  samp.pxm[1] <- in_val
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
      samp.pxm[i] <- propval
      accept <- accept + 1
    }
    else
    {
      samp.pxm[i] <- in_val
    }
    in_val <- samp.pxm[i]
  }
  acceptance <- accept / iter
  print(paste("Acceptance rate is = ", acceptance))
  return(samp.pxm)
}

mybarker <- function(in_val, iter, lambda, delta)
{
  samp.bark <- numeric(length = iter)
  samp.bark[1] <- in_val
  accept <- 0
  for (i in 2:iter)
  {
    propval <- bark.prop(in_val, delta, lambda)
    targ_val.pr <- - prox_arg(propval, prox_func(propval, lambda), lambda)
    targ_val.in <- - prox_arg(in_val, prox_func(in_val, lambda), lambda)
    mh.ratio <- targ_val.pr + log(bark.dens(propval, in_val, delta, lambda)) - 
            targ_val.in - log(bark.dens(in_val, propval, delta, lambda))
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
  acceptance <- accept / iter
  print(paste("Acceptance rate is = ", acceptance))
  return(samp.bark)
}


iter <- 1e5
in_val <- 1
delta <- 1
lambda.vec <- c(1, 0.5, 0.005, 0.00001)

mala.l1 <- mymala(in_val, iter, lambda.vec[1], delta)
mala.l2 <- mymala(in_val, iter, lambda.vec[2], delta)
mala.l3 <- mymala(in_val, iter, lambda.vec[3], delta)
mala.l4 <- mymala(in_val, iter, lambda.vec[4], delta)
mala <- list(mala.l1, mala.l2, mala.l3, mala.l4)


barker.l1 <- mybarker(in_val, iter, lambda.vec[1], delta)
barker.l2 <- mybarker(in_val, iter, lambda.vec[2], delta)
barker.l3 <- mybarker(in_val, iter, lambda.vec[3], delta)
barker.l4 <- mybarker(in_val, iter, lambda.vec[4], delta)
barker <- list(barker.l1, barker.l2, barker.l3, barker.l4)

sam <- seq(-4, 4, length = iter)
den <- exp(-(sam^4)) / (2 * gamma(5/4))

# Density plots

par(mfrow = c(2,2))

for (k in 1:4)
{
  plot(sam, den, type = 'l', xlab = "values",
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

