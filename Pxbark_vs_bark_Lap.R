## Code for PxBarker and Barker algorithms with standard Laplace target

set.seed(123)
library(mcmcse)
library(extraDistr)

target_val <- function(x)
{
  u <- dlaplace(x, 0, 1, log = TRUE)
  return(u)
}

prox_arg <- function(x, y, mu)
{
  z <- abs(y) + ((x-y)^2)/(2*mu)   # proximal function
  return(z)
}

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

tru_log_gradpi <- function(p)
{
  grad <- - p / abs(p)
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
    power_rat <- abs(propval) - abs(in_val)
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
delta <- 2
lambda.vec <- c(1, 0.00001)

# PxBarker samples

pxbarker.l1 <- px.barker(in_val, iter, lambda.vec[1], delta)
pxbarker.l2 <- px.barker(in_val, iter, lambda.vec[2], delta)
pxbarker <- list(pxbarker.l1, pxbarker.l2)

# Barker samples
barker.l1 <- barker_alg(delta, iter, in_val)
barker.l2 <- barker_alg(delta, iter, in_val)
barker <- list(barker.l1, barker.l2)

# Density plots
pdf("plotlap_bark.pdf", height = 6, width = 12)
par(mfrow = c(1,2))

for (k in 1:2)
{
  plot(density(sample <- rlaplace(1e5, 0, 1)), xlab = "values", ylab = "density",
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
  legend("topright", c("PxBarker", "Barker"), lty = 1,
         col = c("blue", "red"), cex = 1, bty = "n")
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

