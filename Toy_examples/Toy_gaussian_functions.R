
############################################################################################
############  Code for Px-MALA and Px-Barker algorithms for pi(x) = k.exp(-x^2) ############
############################################################################################

library(SimTools)

# function calculates the inside of the proximal function
prox_arg <- function(x, y, mu)     # x is the current value
{
  z_x <-  y^2 + ((x-y)^2)/(2*mu)  
  return(z_x)
}

# function calculates the value of the proximal function
prox_func <- function(val, mu)
{
  value <- val/(2*mu + 1) 
  return(value)
}

grad_logpi <- function(val, mu)
{
  gradval <- - (val - prox_func(val, mu)) / mu 
  return(gradval)
}

mymala <- function(curr_val, iter, lambda, delta)
{
  samp.mym <- numeric(length = iter)
  wts_is <- numeric(length = iter)
  
  samp.mym[1] <- curr_val
  proxval_curr <- prox_func(curr_val, lambda)
  targ_val.curr <- - prox_arg(curr_val, proxval_curr, lambda)
  wts_is[1] <-  prox_arg(curr_val, proxval_curr, lambda) - curr_val^2 
  
  accept <- 0
  for (i in 2:iter) 
  {
    prop.mean <- curr_val + (delta / 2)*grad_logpi(curr_val, lambda)
    propval <- prop.mean + sqrt(delta)*rnorm(1)   # proposal step
    
    proxval_next <- prox_func(propval, lambda)
    targ_val.next <- - prox_arg(propval, proxval_next, lambda)
    
    mh.ratio <- targ_val.next + dnorm(curr_val, propval + (delta / 2)*grad_logpi(propval, lambda), 
                                      sqrt(delta), log = TRUE) - targ_val.curr - 
      dnorm(propval, prop.mean, sqrt(delta), log = TRUE)
    if(log(runif(1)) <= mh.ratio)
    {
      samp.mym[i] <- propval
      accept <- accept + 1
      proxval_curr <- proxval_next
      targ_val.curr <- targ_val.next
      
      # weights
      wts_is[i] <- - targ_val.curr - propval^2
    }
    else
    {
      samp.mym[i] <- curr_val
      wts_is[i] <-  - targ_val.curr - curr_val^2 
    }
    curr_val <- samp.mym[i]
  }
  acceptance <- accept / iter
  print(paste("Acceptance rate is = ", acceptance))
  output <- list(samp.mym, wts_is)
  return(output)
}

pxmala <- function(curr_val, iter, lambda, delta)
{
  samp.pxm <- numeric(length = iter)
  
  samp.pxm[1] <- curr_val
  targ_val.curr <- - curr_val^2
  
  accept <- 0
  for (i in 2:iter) 
  {
    prop.mean <- curr_val + (delta / 2)*grad_logpi(curr_val, lambda)
    propval <- prop.mean + sqrt(delta)*rnorm(1)   # proposal step
    
    targ_val.next <- - propval^2
    
    mh.ratio <- targ_val.next + dnorm(curr_val, propval + (delta / 2)*grad_logpi(propval, lambda), 
                                      sqrt(delta), log = TRUE) - targ_val.curr - 
      dnorm(propval, prop.mean, sqrt(delta), log = TRUE)
    if(log(runif(1)) <= mh.ratio)
    {
      samp.pxm[i] <- propval
      accept <- accept + 1
      targ_val.curr <- targ_val.next
    }
    else
    {
      samp.pxm[i] <- curr_val
    }
    curr_val <- samp.pxm[i]
  }
  acceptance <- accept / iter
  print(paste("Acceptance rate is = ", acceptance))
  return(samp.pxm)
}

opt_lambda_func <- function(log_weights)
{
  true_wts <- exp(log_weights)
  ratio <- (mean(true_wts)^2)/mean(true_wts^2)
  return(ratio)
}

# to evaluate asymptotic covariance matrix
asymp_covmat_fn <- function(chain, weights)
{
  wts_mean <- mean(weights)
  num <- chain*weights
  is_est <- sum(num) / sum(weights)
  input_mat <- cbind(num, weights)  # input samples for mcse
  Sigma_mat <- mcse.multi(input_mat)$cov  # estimated covariance matrix of the tuple
  kappa_eta_mat <- cbind(1/wts_mean, -is_est/wts_mean) # derivative of kappa matrix
  asymp_covmat <- (kappa_eta_mat %*% Sigma_mat) %*% t(kappa_eta_mat) # IS asymptotic variance
  return(asymp_covmat)
}

###########################################################################
##################  Barker and other old visualisations ###################
###########################################################################
# bark.prop <- function(val, delta, lambda)
# {
#   y <- rnorm(1, 0, sqrt(delta))
#   prob <- 1 / (1 + exp( -y*grad_logpi(val, lambda)))
#   ifelse(runif(1) <= prob, prop <- val + y, prop <- val - y)
#   return(prop)
# }
# 
# bark.dens <- function(curr_val, propval, delta, lambda)
# {
#   numer <- 2 * dnorm(propval, curr_val, sqrt(delta))
#   denom <- 1 + exp( - (propval - curr_val) * grad_logpi(curr_val, lambda))
#   value <- numer / denom
#   return(value)
# }


# delta.vec <- c(10, 50, 100)

# 
# mala.del1 <- mymala(curr_val, iter, lambda.vec[2], delta.vec[1])
# mala.del2 <- mymala(curr_val, iter, lambda.vec[2], delta.vec[2])
# mala.del3 <- mymala(curr_val, iter, lambda.vec[2], delta.vec[3])
# mala.del <- list(mala.del1, mala.del2, mala.del3)
# 
# pdf("plotexp_lam_MY.pdf", height = 6, width = 15)
#
# # Density plots for varying delta
# 
# par(mfrow = c(1,3))
# 
# for (k in 1:3)
# {
#   plot(sam, den, type = 'l', xlab = "values", ylab = "density",
#        main = bquote(h == .(delta.vec[k])), col = "black")
#   lines(density(as.numeric(unlist(mala.del[k]))) ,col = "blue")
#   lines(density(as.numeric(unlist(barker.del[k]))) ,col = "red")
#   legend("topright", c("Truth", "MYMALA", "MYBarker"), lty = 1,
#          col = c("black", "blue", "red"), cex = 1, bty = "n")
# }
# 
# # acf plots
# 
# par(mfrow = c(1,3))
# 
# for (k in 1:3)
# {
#   acf.mala <- acf(as.numeric(unlist(mala.del[k])), plot = FALSE)$acf
#   acf.barker <- acf(as.numeric(unlist(barker.del[k])), plot = FALSE)$acf
#   plot(1:length(acf.mala), acf.mala, col = "blue", type = 'l',
#        xlab = "Lag", ylab = "Autocorrelation", ylim = c(- .05, 1))
#   lines(1:length(acf.mala), acf.barker, col = "red", type = 'l')
#   legend("topright", c("Truth", "MYMALA", "MYBarker"), lty = 1,
#          col = c("black", "blue", "red"), cex = 1, bty = "n")
# }
# 
# # Trace plots for varying delta
# 
# par(mfrow = c(1,3))
# 
# for (k in 1:3)
# {
#   plot.ts(as.numeric(unlist(barker.del[k])[1:1e3]), col = "red", ylab = "sample values",
#           ylim = c(-5, 5))
#   lines(as.numeric(unlist(mala.del[k])[1:1e3]), col = "blue")
#   legend("topright", c("MYMALA", "MYBarker"), lty = 1,
#          col = c("blue", "red"), cex = 1, bty = "n")
# }
# dev.off()
# mybarker <- function(curr_val, iter, lambda, delta)
# {
#   samp.bark <- numeric(length = iter)
#   samp.bark[1] <- curr_val
#   accept <- 0
#   for (i in 2:iter)
#   {
#     propval <- bark.prop(curr_val, delta, lambda)
#     targ_val.next <- - prox_arg(propval, prox_func(propval, lambda), lambda)
#     targ_val.curr <- - prox_arg(curr_val, prox_func(curr_val, lambda), lambda)
#     mh.ratio <- targ_val.next + log(bark.dens(propval, curr_val, delta, lambda)) - 
#             targ_val.curr - log(bark.dens(curr_val, propval, delta, lambda))
#     if(log(runif(1)) <= mh.ratio)
#     {
#       samp.bark[i] <- propval
#       accept <- accept + 1
#     }
#     else
#     {
#       samp.bark[i] <- curr_val
#     }
#     curr_val <- samp.bark[i]
#   }
#   acceptance <- accept / iter
#   print(paste("Acceptance rate is = ", acceptance))
#   return(samp.bark)
# }

# # MYBarker samples
# 
# barker.l1 <- mybarker(curr_val, iter, lambda.vec[1], delta)
# barker.l2 <- mybarker(curr_val, iter, lambda.vec[2], delta)
# barker <- list(barker.l1, barker.l2)


# 
# barker.del1 <- mybarker(curr_val, iter, lambda.vec[2], delta.vec[1])
# barker.del2 <- mybarker(curr_val, iter, lambda.vec[2], delta.vec[2])
# barker.del3 <- mybarker(curr_val, iter, lambda.vec[2], delta.vec[3])
# barker.del <- list(barker.del1, barker.del2, barker.del3)

