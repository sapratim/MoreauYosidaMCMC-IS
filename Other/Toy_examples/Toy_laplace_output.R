#####################################################################
################# Exponential output visualisation #####################
#####################################################################

library(mcmcse)
source("Toy_laplace_functions.R")
iter <- 1e4
curr_val <- 1

##########  opt lambda visualisation  ######################
lamb <- seq(0.1, 8, length = 1e3)
delta_is <- seq(3, 20, length = length(lamb))
delta_px <- seq(2.8, 1.5, length = length(lamb))

asymp_var_is <- numeric(length = length(lamb))
asymp_var_px <- numeric(length = length(lamb))
rel_eff <- numeric(length = length(lamb))
ratio_vals <- numeric(length = length(lamb))

for(i in 1:length(lamb))
{
  mym.output <- mymala(curr_val, iter, lamb[i], delta_is[i])
  pxm.output <- pxmala(curr_val, iter, lamb[i], delta_px[i])
  asymp_var_is[i] <- asymp_covmat_fn(mym.output[[1]],exp(mym.output[[1]]))
  asymp_var_px[i] <- mcse.multi(pxm.output)$cov
  rel_eff[i] <- asymp_var_px[i]/asymp_var_is[i]
  ratio_vals[i] <- opt_lambda_func(mym.output[[2]])
}

pdf(file = "ess_plot_laplace.pdf", height = 8, width = 10)
par(mfrow = c(2,2))

plot(lamb, ratio_vals, type = 'l', xlab = "lambda", ylab = "imp_ess", 
     main = "ESS for different lambda")
lines(abline(h=c(0.4,0.6), col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))
lines(abline(v=c(max(lamb[which(ratio_vals >= 0.4)]), min(lamb[which(ratio_vals <= 0.6)])), 
             col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb, asymp_var_is, type = 'l', xlab = "lambda", ylab = "asymptotic variance", 
     main = "Asymptotic variance for different lambda")
lines(abline(v=c(max(lamb[which(ratio_vals >= 0.4)]), min(lamb[which(ratio_vals <= 0.6)])), 
             col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb, rel_eff, type = 'l', xlab = "lambda", ylab = "relative efficiency", 
     main = "Relative efficiency for different lambda")
lines(abline(v=c(max(lamb[which(ratio_vals >= 0.4)]), min(lamb[which(ratio_vals <= 0.6)])), 
             col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))
dev.off()

############################################################

lambda.val <- 1
delta_is <- 4
delta_px <- 0.7

#  MYMALA samples 
mala_is <- mymala(in_val, iter, lambda.vec, delta <- 2.8)
mala_px <- pxmala(in_val, iter, lambda.vec, delta_px <- 2.8)
mala_isdens <- mymala(curr_val, iter, lambda.val <- 1e-4, delta_is <- 1.26)  ##  for density plots

# Density plots
pdf("density_trial.pdf", height = 6, width = 12)
par(mfrow=c(1,2))

plot(density(sample <- rlaplace(1e5, 0, 1)), xlab = "values",
     ylab = "density", main = bquote(lambda == .(lambda.vec[k])), col = "black")
lines(density(as.numeric(unlist(mala_is[[1]]))) ,col = "blue")
legend("topright", c("Truth", "MYMALA"), lty = 1,
       col = c("black", "blue"), cex = 1, bty = "n")

plot(density(sample <- rlaplace(1e5, 0, 1)), xlab = "values",
     ylab = "density", main = bquote(lambda == .(lambda.vec[k])), col = "black")
lines(density(as.numeric(unlist(mala_px))) ,col = "blue")
legend("topright", c("Truth", "MYMALA"), lty = 1,
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
