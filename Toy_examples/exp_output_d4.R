
source("Toy_exp_functions.R")
load("output_d4.Rdata")
curr_val <- 1

reps <- 50
cutoff <- c(0.4,0.8)
lamb_d4 <- seq(1e-4, 25, length = reps)
load("output_d4.Rdata")

ess_d4 <- sapply(output_exp_d4, function(l) l[[4]])
pilam_var_d4 <- sapply(output_exp_d4, function(l) l[[2]])
is_var_d4 <- sapply(output_exp_d4, function(l) l[[1]])

opt_win <- c(max(lamb_d4[which(ess_d4 >= cutoff[1])]), min(lamb_d4[which(ess_d4 <= cutoff[2])]))

pdf(file = "exp_plot_beta4.pdf", height = 4, width = 12)
par(mfrow = c(1,3))

plot(lamb_d4, ess_d4, type = 'l', xlab = expression(lambda), ylab = expression(n[e]/n))
# lines(abline(h=cutoff, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))
lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb_d4, pilam_var_d4, type = 'l', xlab = expression(lambda), 
     ylab = expression(paste("MCMC ESS under ",pi^lambda))) 
# lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb_d4, is_var_d4, type = 'l', xlab = expression(lambda), 
     ylab ="Importance Sampling Asymptotic variance")
lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))
dev.off()