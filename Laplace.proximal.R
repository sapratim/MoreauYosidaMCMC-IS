
 ## Code for Px-MALA and Px-Barker algorithms with standard Laplace target

library(extraDistr)

target_val <- function(x)
{
  u <- dlaplace(x, 0, 1)
  return(u)
}

proxfunc <- function(x, y, mu)
{
  z <- abs(y) + ((x-y)^2)/(2*mu)   # proximal function
  return(z)
}

px.mala <- function(in_val, iter, lambda, delta)
  {
  frac <- delta / 2*lambda
  noise <- rnorm(iter, 0, delta)
  samp.pxm <- numeric(length = iter)
  samp.pxm[1] <- in_val
  
  for (i in 2:iter) 
  {
    ifelse(in_val + lambda < 0, dummy.in <- in_val + lambda, dummy.in <- in_val - lambda)
    val0.in <- proxfunc(in_val, 0, lambda)
    valdummy.in <- proxfunc(in_val, dummy.in, lambda)
    ifelse(val0.in < valdummy.in, prox.in <- 0, prox.in <- dummy.in)
    propval <- (1 - frac)*in_val + frac*prox.in + sqrt(delta)*noise[i]   # proposal step
    
    ifelse(propval + lambda < 0, dummy.prop <- propval + lambda, dummy.prop <- propval - lambda)
    val0.prop <- proxfunc(propval, 0, lambda)
    valdummy.prop <- proxfunc(propval, dummy.prop, lambda)
    ifelse(val0.prop < valdummy.prop, prox.prop <- 0, prox.prop <- dummy.prop)
    mh.ratio <- (target_val(propval)*dnorm(in_val, (1 - frac)*propval + frac*prox.prop, delta)) /
           (target_val(in_val)*dnorm(propval, (1 - frac)*in_val + frac*prox.in, delta))
      if(runif(1) <= mh.ratio)
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
in_val <- 5
lambda.vec <- c(1, 0.1, 0.01, 0.001)
delta.vec <- rep(1, 4)

plot(density(sample <- rlaplace(1e5, 0, 1)), xlab = "values",
     ylab = "density", main = "Proximal MALA for standard Laplace", col = "black")
lines(density(samp.l1 <- px.mala(in_val, iter, lambda.vec[1], delta.vec[1])) ,col = "blue")
lines(density(samp.l2 <- px.mala(in_val, iter, lambda.vec[2], delta.vec[2])) ,col = "green")
lines(density(samp.l3 <- px.mala(in_val, iter, lambda.vec[3], delta.vec[3])) ,col = "orange")
lines(density(samp.l4 <- px.mala(in_val, iter, lambda.vec[4], delta.vec[4])) ,col = "red")
legend("topright", c("Truth", "lambda = 1", "lambda = 0.1", "lambda = 0.01", "lambda = 0.001"),
   col = c("black", "blue", "green", "orange", "red"), fill = c("black", "blue", "green", "orange", "red"))

