
############################################################################################
############  Code for Px-MALA and Px-Barker algorithms for pi(x) = k.exp(-x^4) ############
############################################################################################

library(SimTools)

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
  wts_is[1] <-  prox_arg(curr_val, proxval_curr, lambda) - curr_val^4 
  
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
      wts_is[i] <- prox_arg(propval, proxval_next, lambda) - propval^4
    }
    else
    {
      samp.mym[i] <- curr_val
      wts_is[i] <-  prox_arg(curr_val, proxval_curr, lambda) - curr_val^4 
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
  targ_val.curr <- - curr_val^4

  accept <- 0
  for (i in 2:iter) 
  {
    prop.mean <- curr_val + (delta / 2)*grad_logpi(curr_val, lambda)
    propval <- prop.mean + sqrt(delta)*rnorm(1)   # proposal step
    
    targ_val.next <- - propval^4
    
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

opt_lambda_func <- function(curr_val, iter, lambda, delta)
{
  mym.output <- mymala(curr_val, iter, lambda, delta)
  true_wts <- exp(mym.output[[2]])
  mean_wts <- mean(true_wts)^2
  mean_sq_wts <- mean(true_wts^2)
  ratio <- (mean_wts/mean_sq_wts)
  return(ratio)
}

################################################################################
##############################   Visualisations ################################
################################################################################


iter <- 1e5
curr_val <- 1
lambda.val <- 1
delta_is <- 4
delta_px <- 0.7

##########  opt lambda visualisation  ######################
lamb <- seq(0.0001, 100, length = 200)
ratio_vals <- numeric(length = length(lamb))
for(i in 1:length(lamb))
{
  ratio_vals[i] <- opt_lambda_func(curr_val <- 1, iter <- 1e5, lambda <- lamb[i], delta <- 4)
}
pdf(file = "ess_plot.pdf", height = 8, width = 10)
plot(lamb, ratio_vals, type = 'l', xlab = "lambda", ylab = "imp_ess", 
     main = "ESS for different lambda")
lines(abline(h=c(0.4,0.6), col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))
dev.off()
############################################################

#  MYMALA samples 
mala.is <- mymala(curr_val, iter, lambda.val, delta_is)
mala.px <- pxmala(curr_val, iter, lambda.val, delta_px)
mala_isdens <- mymala(curr_val, iter, lambda.val <- 1e-4, delta_is <- 0.75)  ##  for density plots

# Actual density shape
sam <- seq(-4, 4, length = iter)
den <- exp(-(sam^4)) / (2 * gamma(5/4))

# Density plots
pdf("density_trial.pdf", height = 6, width = 12)
par(mfrow=c(1,2))

plot(sam, den, type = 'l', xlab = "values", ylab = "density",
        col = "black")
lines(density(as.numeric(unlist(mala_isdens[1]))) ,col = "blue")
legend("topright", c("Truth", "MYMALA"), lty = 1,
       col = c("black", "blue"), cex = 1, bty = "n")
  
plot(sam, den, type = 'l', xlab = "values", ylab = "density",
       col = "black")
lines(density(mala.px) ,col = "blue")
legend("topright", c("Truth", "PxMALA"), lty = 1,
       col = c("black", "blue"), cex = 1, bty = "n")
dev.off()

# acf plots
pdf("acf_trial.pdf", height = 8, width = 10)

  acf.ism <- acf(as.numeric(unlist(mala.is[1])), plot = FALSE)$acf
  acf.pxm <- acf(mala.px, plot = FALSE)$acf
  plot(1:length(acf.ism), acf.ism, col = "blue", type = 'l',
       xlab = "Lag", ylab = "Autocorrelation", ylim = c(-0.25, 1))
  lines(1:length(acf.pxm), acf.pxm, col = "red", type = 'l')
  legend("topright", c("MYMALA", "PxMALA"), lty = 1,
         col = c("blue", "red"), cex = 1, bty = "n")
dev.off()

# Trace plots
pdf("trace_trial.pdf", height = 8, width = 10)
traceplot(list(mala.is[[1]], mala.px), col = c("blue", "red"))
dev.off()

#################  Boxplots of weights ################# 
wt_mat <- matrix(0, nrow = 1e4, ncol = length(lamb))

for (j in 1:length(lamb)) {
  wt_mat[,j] <- mymala(curr_val, iter<-1e4, lamb[j], delta <- 2)[[2]]
}
pdf(file = "boxplots_wts_trial.pdf", width = 12, height = 8)
boxplot(exp(wt_mat), use.cols = TRUE, xlab = "Coordinate", 
        ylab = "Average relative efficiency")
dev.off()

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
