
###########  PxBarker vs Barker for exp(-x^4)

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
  z_x <-  y^4 + ((x-y)^2)/2*mu  
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

log_gradpi <- function(val, mu)    #### gradient of the MY envelope
{
  gradval <- - (val - prox_func(val, mu)) / mu 
  return(gradval)
}

tru_log_gradpi <- function(p)
{
  grad <- - 4*(p^3)
  return(grad)
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
  acceptance <- accept/iter
  print(paste("Acceptance rate is = ", acceptance))
  return(samp.bark)
}


barker_alg <- function(delta, iter, in_val)
{
  samp.bark <- numeric(length = iter)
  samp.bark[1] <- in_val
  accept <- 0
  # iterations using Barker's algorithm
  err <- rnorm(iter, 0, sqrt(delta))
  U <- runif(iter)
  V <- runif(iter)
  for (i in 2:iter) 
    {
    z <- err[i]
    prob <- 1 / (1 + exp(- z * tru_log_gradpi(in_val)))
    propval <- ifelse(U[i] < prob, in_val + z, in_val - z)
    ratio_prop = (1+exp((in_val - propval)*tru_log_gradpi(in_val))) /
                  (1+exp((propval - in_val)*tru_log_gradpi(propval)))
    power_rat <- propval^4 - in_val^4
    acc_prob <- min(1, exp(- power_rat)*ratio_prop)
    if(V[i] < acc_prob)
    {
      samp.bark[i] = propval
      accept = accept + 1
    }
    else
    { 
      samp.bark[i] = in_val
    }
    in_val = samp.bark[i]
  }
  acceptance <- accept/iter
  print(paste("Acceptance rate is = ", acceptance))
  return(samp.bark)
}
  
iter <- 1e5
in_val <- 1
delta <- 10
lambda.vec <- c(1, 0.00001)

# PxBarker samples

pxbarker.l1 <- px.barker(in_val, iter, lambda.vec[1], delta)
pxbarker.l2 <- px.barker(in_val, iter, lambda.vec[2], delta)
pxbarker <- list(pxbarker.l1, pxbarker.l2)

# Barker samples
barker.l1 <- barker_alg(delta, iter, in_val)
barker.l2 <- barker_alg(delta, iter, in_val)
barker <- list(barker.l1, barker.l2)

# Actual density shape
sam <- seq(-4, 4, length = iter)
den <- exp(-(sam^4)) / (2 * gamma(5/4))

# Density plots
pdf("plotexp_bark.pdf", height = 6, width = 12)
par(mfrow = c(1,2))

for (k in 1:2)
{
  plot(sam, den, type = 'l', xlab = "values", ylab = "density",
       main = bquote(lambda == .(lambda.vec[k])), col = "black")
  lines(density(as.numeric(unlist(pxbarker[k]))) ,col = "blue")
  lines(density(as.numeric(unlist(barker[k]))) ,col = "red")
  legend("topright", c("Truth", "PxBarker", "Barker"), lty = 1,
         col = c("black", "blue", "red"), cex = 1, bty = "n")
}

# acf plots

par(mfrow = c(1,2))

for (k in 1:2)
{
  acf.pxbarker <- acf(as.numeric(unlist(pxbarker[k])), plot = FALSE)$acf
  acf.barker <- acf(as.numeric(unlist(barker[k])), plot = FALSE)$acf
  plot(1:length(acf.pxbarker), acf.pxbarker, col = "blue", type = 'l',
       xlab = "Lag", ylab = "Autocorrelation", ylim = c(-0.2, 1))
  lines(1:length(acf.pxbarker), acf.barker, col = "red", type = 'l')
  legend("topright", c("Truth", "PxBarker", "Barker"), lty = 1,
         col = c("black", "blue", "red"), cex = 1, bty = "n")
}

# Trace plots

par(mfrow = c(1,2))

for (k in 1:2)
{
  plot.ts(as.numeric(unlist(barker[k])[1:1e3]), col = "red", ylab = "sample values",
          ylim = c(-5, 5))
  lines(as.numeric(unlist(pxbarker[k])[1:1e3]), col = "blue")
  legend("topright", c("PxBarker", "Barker"), lty =1,
         col = c("blue", "red"), cex = 1, bty = "n")
}
dev.off()
