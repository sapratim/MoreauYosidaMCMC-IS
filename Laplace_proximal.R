
## Code for Px-MALA and Px-Barker algorithms with standard Laplace target

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

px.mala <- function(in_val, iter, lambda, delta)
  {
  samp.pxm <- numeric(length = iter)
  samp.pxm[1] <- in_val
  
  for (i in 2:iter) 
  {
    propval <- rnorm(1, in_val + (delta / 2)*log_gradpi(in_val, lambda), sqrt(delta))   # proposal step
    mh.ratio <- target_val(propval) + dnorm(in_val, propval + (delta / 2)*log_gradpi(propval, lambda), delta, log = TRUE) - 
           target_val(in_val) - dnorm(propval, in_val + (delta / 2)*log_gradpi(in_val, lambda), delta, log = TRUE)
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

iter <- 1e5
in_val <- 1
delta.vec <- rep(1, 4)
lambda.vec <- c(1, 0.5, 0.005, 0.00001)

plot(density(sample <- rlaplace(1e5, 0, 1)), xlab = "values",
     ylab = "density", main = "Proximal MALA for standard Laplace", col = "black")
lines(density(samp.l1 <- px.mala(in_val, iter, lambda.vec[1], delta.vec[1])) ,col = "blue")
lines(density(samp.l2 <- px.mala(in_val, iter, lambda.vec[2], delta.vec[2])) ,col = "green")
lines(density(samp.l3 <- px.mala(in_val, iter, lambda.vec[3], delta.vec[3])) ,col = "orange")
lines(density(samp.l4 <- px.mala(in_val, iter, lambda.vec[4], delta.vec[4])) ,col = "red")
legend("topright", c("Truth", "lambda = 1", "lambda = 0.1", "lambda = 0.01", "lambda = 0.001"),
   col = c("black", "blue", "green", "orange", "red"), fill = c("black", "blue", "green", "orange", "red"))

# acf plots
par(mfrow = c(2,2))
acf(samp.l1)
acf(samp.l2)
acf(samp.l3)
acf(samp.l4)

# Trace plots
par(mfrow = c(2,2))
plot.ts(samp.l1)
plot.ts(samp.l2)
plot.ts(samp.l3)
plot.ts(samp.l4)
