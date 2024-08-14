

library(mcmcse)
pdf(file = "ani_plot_d5.pdf", height = 8, width = 10)
par(mfrow = c(2,2))

plot(lamb_d5, ess_ism_d5, type = 'l', xlab = "lambda", ylab = "imp_ess", 
     main = "ESS for different lambda")
lines(abline(h=c(0.4,0.6), col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))
lines(abline(v=c(max(lamb_d5[which(ess_ism_d5 >= 0.4)]), min(lamb_d5[which(ess_ism_d5 <= 0.6)])), 
             col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb_d5, asymp_var_pilam, type = 'l', xlab = "lambda", ylab = "asymptotic variance", 
     main = paste("Asymptotic variance  for pi^{lambda}"))
lines(abline(v=c(max(lamb_d5[which(ess_ism_d5 >= 0.4)]), min(lamb_d5[which(ess_ism_d5 <= 0.6)])), 
             col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb_d5, asymp_margvar_ism_d5, type = 'l', xlab = "lambda", ylab = "asymptotic variance", 
     main = "Asymptotic variance for pi based on IS estimator")
lines(abline(v=c(max(lamb_d5[which(ess_ism_d5 >= 0.4)]), min(lamb_d5[which(ess_ism_d5 <= 0.6)])), 
             col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))
dev.off()






pdf(file = "ani_plot_d10.pdf", height = 8, width = 10)
par(mfrow = c(2,2))

plot(lamb, ess_ism_d10, type = 'l', xlab = "lambda", ylab = "imp_ess", 
     main = "ESS for different lambda")
lines(abline(h=c(0.4,0.6), col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))
lines(abline(v=c(max(lamb[which(ess_ism_d10 >= 0.4)]), min(lamb[which(ess_ism_d10 <= 0.6)])), 
             col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb, asymp_margvar_ism_d10, type = 'l', xlab = "lambda", ylab = "asymptotic variance", 
     main = "Asymptotic variance for different lambda")
lines(abline(v=c(max(lamb[which(ess_ism_d10 >= 0.4)]), min(lamb[which(ess_ism_d10 <= 0.6)])), 
             col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb, avg_rel_eff_10, type = 'l', xlab = "lambda", ylab = "relative efficiency", 
     main = "Relative efficiency for different lambda")
lines(abline(v=c(max(lamb[which(ess_ism_d10 >= 0.4)]), min(lamb[which(ess_ism_d10 <= 0.6)])), 
             col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))
dev.off()






pdf(file = "ani_plot_d50.pdf", height = 8, width = 10)
par(mfrow = c(2,2))

plot(lamb_d50, ess_ism_d50, type = 'l', xlab = "lambda", ylab = "imp_ess", 
     main = "ESS for different lambda")
lines(abline(h=c(0.4,0.6), col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))
lines(abline(v=c(max(lamb_d50[which(ess_ism_d50 >= 0.4)]), min(lamb_d50[which(ess_ism_d50 <= 0.6)])), 
             col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb_d50, asymp_var_pilam_d50, type = 'l', xlab = "lambda", ylab = "asymptotic variance", 
     main = "Asymptotic variance for pi^{lambda}")
lines(abline(v=c(max(lamb_d50[which(ess_ism_d50 >= 0.4)]), min(lamb_d50[which(ess_ism_d50 <= 0.6)])), 
             col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb_d50, asymp_margvar_ism_d50, type = 'l', xlab = "lambda", ylab = "asymptotic variance", 
     main = "Asymptotic variance for pi based on IS estimator")
lines(abline(v=c(max(lamb_d50[which(ess_ism_d50 >= 0.4)]), min(lamb_d50[which(ess_ism_d50 <= 0.6)])), 
             col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))
dev.off()

