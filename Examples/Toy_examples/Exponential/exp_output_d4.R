
source("exp_d4_functions.R")
load("output_d1.Rdata")
load("output_d5.Rdata")
load("output_d10.Rdata")
load("output_d20.Rdata")

reps <- 50
cutoff <- c(0.4,0.8)
lamb_d1 <- seq(1e-4, 5, length = reps)
lamb_d5 <- seq(1e-4, .6, length = reps)
lamb_d10 <- seq(1e-3, 0.5, length = reps)
lamb_d20 <- seq(1e-3, 0.2, length = reps)


##########################  For d = 1  ##########################
load("output_d1.Rdata")

ess_d1 <- sapply(output_exp_d1, function(l) l[[4]])
pilam_var_d1 <- sapply(output_exp_d1, function(l) l[[2]])
is_var_d1 <- sapply(output_exp_d1, function(l) l[[1]])

opt_win <- c(max(lamb_d1[which(ess_d1 >= cutoff[1])]), min(lamb_d1[which(ess_d1 <= cutoff[2])]))

pdf(file = "plots/exp_plot_d1.pdf", height = 4, width = 12)
par(mfrow = c(1,3))

plot(lamb_d1, ess_d1, type = 'l', xlab = expression(lambda), ylab = expression(n[e]/n))
# lines(abline(h=cutoff, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))
lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb_d1, pilam_var_d1, type = 'l', xlab = "lambda", 
     ylab = expression(paste("MCMC ESS under ",pi^lambda))) 
# lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb_d1, is_var_d1, type = 'l', xlab = "lambda",
     ylab ="Importance Sampling Asymptotic variance")
lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))
dev.off()

##########################  For d = 5  ##########################
load("output_d5.Rdata")

ess_d5 <- sapply(output_exp_d5, function(l) l[[4]])
pilam_var_d5 <- colMeans(sapply(output_exp_d5, function(l) l[[2]]))
is_var_d5 <- colMeans(sapply(output_exp_d5, function(l) l[[1]]))

opt_win <- c(max(lamb_d5[which(ess_d5 >= cutoff[1])]), min(lamb_d5[which(ess_d5 <= cutoff[2])]))


pdf(file = "plots/exp_plot_d5.pdf", height = 4, width = 12)
par(mfrow = c(1,3))

plot(lamb_d5, ess_d5, type = 'l', xlab = expression(lambda), ylab = expression(n[e]/n))
# lines(abline(h=cutoff, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))
lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb_d5, pilam_var_d5, type = 'l', xlab = expression(lambda), 
     ylab = expression(paste("MCMC ESS under ",pi^lambda))) 
# lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb_d5, is_var_d5, type = 'l', xlab = expression(lambda), 
     ylab ="Importance Sampling Asymptotic variance")

lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))
dev.off()


##########################  For d = 10  ##########################
load("output_d10.Rdata")

ess_d10 <- sapply(output_exp_d10, function(l) l[[4]])
pilam_var_d10 <- colMeans(sapply(output_exp_d10, function(l) l[[2]]))
is_var_d10 <- colMeans(sapply(output_exp_d10, function(l) l[[1]]))

opt_win <- c(max(lamb_d10[which(ess_d10 >= cutoff[1])]), min(lamb_d10[which(ess_d10 <= cutoff[2])]))

pdf(file = "plots/exp_plot_d10.pdf", height = 4, width = 12)
par(mfrow = c(1,3))

plot(lamb_d10, ess_d10, type = 'l', xlab = expression(lambda), ylab = expression(n[e]/n))
# lines(abline(h=cutoff, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))
lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb_d10, pilam_var_d10, type = 'l', xlab = expression(lambda), 
     ylab = expression(paste("MCMC ESS under ",pi^lambda))) 
# lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb_d10, is_var_d10, type = 'l', xlab = expression(lambda), 
     ylab ="Importance Sampling Asymptotic variance")
lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))
dev.off()


##########################  For d = 20  ##########################
load("output_d20.Rdata")

ess_d20 <- sapply(output_exp_d20, function(l) l[[4]])
pilam_var_d20 <- colMeans(sapply(output_exp_d20, function(l) l[[2]]))
is_var_d20 <- colMeans(sapply(output_exp_d20, function(l) l[[1]]))

opt_win <- c(max(lamb_d20[which(ess_d20 >= cutoff[1])]), min(lamb_d20[which(ess_d20 <= cutoff[2])]))


pdf(file = "plots/exp_plot_d20.pdf", height = 4, width = 12)
par(mfrow = c(1,3))

plot(lamb_d20, ess_d20, type = 'l', xlab = expression(lambda), ylab = expression(n[e]/n))
# lines(abline(h=cutoff, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))
lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb_d20, pilam_var_d20, type = 'l', xlab = expression(lambda), 
     ylab = expression(paste("MCMC ESS under ",pi^lambda))) 
# lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb_d20, is_var_d20, type = 'l', xlab = expression(lambda), 
     ylab ="Importance Sampling Asymptotic variance")
lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))
dev.off()
