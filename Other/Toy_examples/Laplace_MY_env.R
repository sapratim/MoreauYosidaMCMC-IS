
###############   Function for the MY envelope of Laplace distribution   #############

values <- seq(-3, 3, length = 1e5)

prox_arg <- function(x, y, mu)     # x is the current value
{
  z_x <- abs(y) + (x - y)^2 / (2*mu)    
  return(z_x)
}

prox_func <- function(val, mu)
{
  vec <- c(0, val - mu, val + mu)
  fun.val <- prox_arg(val, vec, mu)
  index <- which.min(fun.val)
  prox <- vec[index]
  return(prox)
}

approx_dens_func <- function(x, lambda)   ########   from Durmus et al.
{
  num <- exp((lambda/2 - abs(x))*as.numeric(I(abs(x) >= lambda)) - 
               ((x^2)/(2*lambda))*as.numeric(I(abs(x) < lambda)))
  den <- 2*(exp(-lambda/2) + (((2*pi)*lambda)^(1/2))*(pnorm(lambda^(1/2)) - .5))
  val <- num/den
  return(val)
}
      
psi_vals_true <- abs(values)
psi_vals_1 <- prox_arg(values, sapply(values, function (x) prox_func(x,mu = 1)), mu = 1)
psi_vals_0.5 <- prox_arg(values, sapply(values, function (x) prox_func(x,mu = 0.5)), mu = 0.5)
psi_vals_0.25 <- prox_arg(values, sapply(values, function (x) prox_func(x,mu = 0.25)), mu = 0.25)

true_dens_vals <- exp(-abs(values))/2   #######  corresponding to lambda = 0
approx_dens_1 <- approx_dens_func(values, lambda = 1)
approx_dens_0.5 <- approx_dens_func(values, lambda = 0.5)
approx_dens_0.25 <- approx_dens_func(values, lambda = 0.25)

pdf(file = "MY_env_lap.pdf", height = 6, width = 12)
par(mfrow = c(1,2))

plot(values, psi_vals_true, type = 'l', xlab = "x", ylab = "Envelope",
            main = expression(paste(psi^lambda)))
lines(values, psi_vals_0.25, type = 'l', col = "red")
lines(values, psi_vals_0.5, type = 'l', col = "blue")
lines(values, psi_vals_1, type = 'l', col = "magenta")
legend("bottomright", c("0", "0.25", "0.5", "1"), lty = 1,
       col = c("black", "red", "blue", "magenta"), cex = 0.6, bty = "n")

plot(values, true_dens_vals, type = 'l', xlab = "x", ylab = "Density",
     main = expression(paste(pi^lambda)))
lines(values, approx_dens_0.25, type = 'l', col = "red")
lines(values, approx_dens_0.5, type = 'l', col = "blue")
lines(values, approx_dens_1, type = 'l', col = "magenta")
legend("topright", c("0", "0.25", "0.5", "1"), lty = 1,
       col = c("black", "red", "blue", "magenta"), cex = 0.6, bty = "n")
dev.off()
