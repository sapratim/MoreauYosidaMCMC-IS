
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

iter <- 1e5
in_val <- 1
delta.vec <- c(1, 0.5, 0.1, 0.05)
lambda.vec <- delta.vec/2

## Density Plots

plot(density(sample <- rlaplace(1e5, 0, 1)), xlab = "values",
     ylab = "density", ylim = c(0, 0.6), main = "MYMALA for standard Laplace", col = "black")
lines(density(samp.l1 <- px.mala(in_val, iter, lambda.vec[1], delta.vec[1])) ,col = "blue")
lines(density(samp.l2 <- px.mala(in_val, iter, lambda.vec[2], delta.vec[2])) ,col = "green")
lines(density(samp.l3 <- px.mala(in_val, iter, lambda.vec[3], delta.vec[3])) ,col = "orange")
lines(density(samp.l4 <- px.mala(in_val, iter, lambda.vec[4], delta.vec[4])) ,col = "red")
legend("topright", c("Truth", "lambda = 1", "lambda = 0.1", "lambda = 0.01", "lambda = 0.001"),
       col = c("black", "blue", "green", "orange", "red"), fill = c("black", "blue", "green", "orange", "red"))

## Trace plots

par(mfrow = c(2,2))
plot.ts(samp.l1)
plot.ts(samp.l2)
plot.ts(samp.l3)
plot.ts(samp.l4)

## ACF Plots

par(mfrow = c(2,2))
acf(samp.l1)
acf(samp.l2)
acf(samp.l3)
acf(samp.l4)

