
## Code for Px-MALA and Px-Barker algorithms with standard Laplace target

set.seed(123)

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
  acceptance <- accept/iter
  print(paste("Acceptance rate is = ", acceptance))
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
  acceptance <- accept/iter
  print(paste("Acceptance rate is = ", acceptance))
  return(samp.bark)
}


iter <- 1e5
in_val <- 1
delta <- 2
lambda.vec <- c(1, 0.00001)

#  PxMALA samples 

mala.l1 <- px.mala(in_val, iter, lambda.vec[1], delta)
mala.l2 <- px.mala(in_val, iter, lambda.vec[2], delta)
mala <- list(mala.l1, mala.l2)

# PxBarker samples

barker.l1 <- px.barker(in_val, iter, lambda.vec[1], delta)
barker.l2 <- px.barker(in_val, iter, lambda.vec[2], delta)
barker <- list(barker.l1, barker.l2)

# Actual density shape
sam <- seq(-4, 4, length = iter)
den <- exp(-(sam^4)) / (2 * gamma(5/4))

# Density plots
pdf("plotexp_del2.pdf", height = 6, width = 12)
par(mfrow = c(1,2))

for (k in 1:2)
{
  plot(sam, den, type = 'l', xlab = "values", ylab = "density",
       main = bquote(lambda == .(lambda.vec[k])), col = "black")
  lines(density(as.numeric(unlist(mala[k]))) ,col = "blue")
  lines(density(as.numeric(unlist(barker[k]))) ,col = "red")
  legend("topright", c("Truth", "PxMALA", "PxBarker"), lty = 1,
         col = c("black", "blue", "red"), cex = 1, bty = "n")
}

# acf plots

par(mfrow = c(1,2))

for (k in 1:2)
{
  acf.mala <- acf(as.numeric(unlist(mala[k])), plot = FALSE)$acf
  acf.barker <- acf(as.numeric(unlist(barker[k])), plot = FALSE)$acf
  plot(1:length(acf.mala), acf.mala, col = "blue", type = 'l',
       xlab = "Lag", ylab = "Autocorrelation", ylim = c(-0.2, 1))
  lines(1:length(acf.mala), acf.barker, col = "red", type = 'l')
  legend("topright", c("Truth", "PxMALA", "PxBarker"), lty = 1,
         col = c("black", "blue", "red"), cex = 1, bty = "n")
}

# Trace plots

par(mfrow = c(1,2))

for (k in 1:2)
{
  plot.ts(as.numeric(unlist(barker[k])[1:1e3]), col = "red", ylab = "sample values",
          ylim = c(-5, 5))
  lines(as.numeric(unlist(mala[k])[1:1e3]), col = "blue")
  legend("topright", c("PxMALA", "PxBarker"), lty =1,
         col = c("blue", "red"), cex = 1, bty = "n")
}
dev.off()



delta.vec <- c(10, 50, 100)

mala.del1 <- px.mala(in_val, iter, lambda.vec[2], delta.vec[1])
mala.del2 <- px.mala(in_val, iter, lambda.vec[2], delta.vec[2])
mala.del3 <- px.mala(in_val, iter, lambda.vec[2], delta.vec[3])
mala.del <- list(mala.del1, mala.del2, mala.del3)

barker.del1 <- px.barker(in_val, iter, lambda.vec[2], delta.vec[1])
barker.del2 <- px.barker(in_val, iter, lambda.vec[2], delta.vec[2])
barker.del3 <- px.barker(in_val, iter, lambda.vec[2], delta.vec[3])
barker.del <- list(barker.del1, barker.del2, barker.del3)

pdf("plotexp_lam.pdf", height = 6, width = 15)

# Density plots for varying delta

par(mfrow = c(1,3))

for (k in 1:3)
{
  plot(sam, den, type = 'l', xlab = "values", ylab = "density", 
           main = bquote(delta == .(delta.vec[k])), col = "black") 
  lines(density(as.numeric(unlist(mala.del[k]))) ,col = "blue")
  lines(density(as.numeric(unlist(barker.del[k]))) ,col = "red")
  legend("topright", c("Truth", "PxMALA", "PxBarker"), lty = 1,
         col = c("black", "blue", "red"), cex = 1, bty = "n")
}

# acf plots

par(mfrow = c(1,3))

for (k in 1:3)
{
  acf.mala <- acf(as.numeric(unlist(mala.del[k])), plot = FALSE)$acf
  acf.barker <- acf(as.numeric(unlist(barker.del[k])), plot = FALSE)$acf
  plot(1:length(acf.mala), acf.mala, col = "blue", type = 'l',
       xlab = "Lag", ylab = "Autocorrelation", ylim = c(- .05, 1))
  lines(1:length(acf.mala), acf.barker, col = "red", type = 'l')
  legend("topright", c("Truth", "PxMALA", "PxBarker"), lty = 1,
         col = c("black", "blue", "red"), cex = 1, bty = "n")
}

# Trace plots for varying delta

par(mfrow = c(1,3))

for (k in 1:3)
{
  plot.ts(as.numeric(unlist(barker.del[k])[1:1e3]), col = "red", ylab = "sample values",
          ylim = c(-5, 5))
  lines(as.numeric(unlist(mala.del[k])[1:1e3]), col = "blue")
  legend("topright", c("PxMALA", "PxBarker"), lty = 1,
         col = c("blue", "red"), cex = 1, bty = "n")
}
dev.off()
