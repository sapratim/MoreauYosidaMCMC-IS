
load("output_ani_laplace.Rdata")

##########################  For d = 5  ##########################

lamb_d5 <- seq(1e-2, 1, length = 1e2)
pdf(file = "ani_plot_d5.pdf", height = 8, width = 10)
par(mfrow = c(2,2))

plot(lamb_d5, output[[1]][[2]], type = 'l', xlab = "lambda", ylab = "imp_ess", 
     main = "ESS for different lambda")
lines(abline(h=c(0.4,0.6), col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))
lines(abline(v=c(max(lamb_d5[which(output[[1]][[2]] >= 0.4)]),
   min(lamb_d5[which(output[[1]][[2]] <= 0.6)])),  col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb_d5, output[[1]][[3]], type = 'l', xlab = "lambda", ylab = "asymptotic variance", 
     main = paste("Asymptotic variance  for pi^{lambda}"))
lines(abline(v=c(max(lamb_d5[which(output[[1]][[2]] >= 0.4)]), 
   min(lamb_d5[which(output[[1]][[2]] <= 0.6)])), col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb_d5, output[[1]][[1]], type = 'l', xlab = "lambda", ylab = "asymptotic variance", 
     main = "Asymptotic variance for pi based on IS estimator")
lines(abline(v=c(max(lamb_d5[which(output[[1]][[2]] >= 0.4)]), 
   min(lamb_d5[which(output[[1]][[2]] <= 0.6)])), col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))
dev.off()


##########################  For d = 10  ##########################


lamb_d10 <- seq(1e-2, 0.5, length = 1e2)
pdf(file = "ani_plot_d10.pdf", height = 8, width = 10)
par(mfrow = c(2,2))

plot(lamb_d10, output[[2]][[2]], type = 'l', xlab = "lambda", ylab = "imp_ess", 
     main = "ESS for different lambda")
lines(abline(h=c(0.4,0.6), col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))
lines(abline(v=c(max(lamb_d10[which(output[[2]][[2]] >= 0.4)]), 
  min(lamb_d10[which(output[[2]][[2]] <= 0.6)])),  col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb_d10, output[[2]][[3]], type = 'l', xlab = "lambda", ylab = "asymptotic variance", 
     main = "Asymptotic variance  for pi^{lambda}")
lines(abline(v=c(max(lamb_d10[which(output[[2]][[2]] >= 0.4)]), 
   min(lamb_d10[which(output[[2]][[2]] <= 0.6)])), col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb_d10, output[[2]][[1]], type = 'l', xlab = "lambda", ylab = "asymptotic variance", 
     main = "Asymptotic variance for pi based on IS estimator")
lines(abline(v=c(max(lamb_d10[which(output[[2]][[2]] >= 0.4)]),
   min(lamb_d10[which(output[[2]][[2]] <= 0.6)])), col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))
dev.off()


##########################  For d = 20  ##########################

lamb_d20 <- seq(5e-3, 0.1, length = 1e2)
pdf(file = "ani_plot_d20.pdf", height = 8, width = 10)
par(mfrow = c(2,2))

plot(lamb_d20, output[[3]][[2]], type = 'l', xlab = "lambda", ylab = "imp_ess", 
     main = "ESS for different lambda")
lines(abline(h=c(0.4,0.6), col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))
lines(abline(v=c(max(lamb_d20[which(output[[3]][[2]] >= 0.4)]), 
   min(lamb_d20[which(output[[3]][[2]] <= 0.6)])),  col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb_d20, output[[3]][[3]], type = 'l', xlab = "lambda", ylab = "asymptotic variance", 
     main = "Asymptotic variance  for pi^{lambda}")
lines(abline(v=c(max(lamb_d20[which(output[[3]][[2]] >= 0.4)]), 
   min(lamb_d20[which(output[[3]][[2]] <= 0.6)])),  col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb_d20, output[[3]][[1]], type = 'l', xlab = "lambda", ylab = "asymptotic variance", 
     main = "Asymptotic variance for pi based on IS estimator")
lines(abline(v=c(max(lamb_d20[which(output[[3]][[2]] >= 0.4)]), 
   min(lamb_d20[which(output[[3]][[2]] <= 0.6)])),  col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))
dev.off()


##########################  For d = 50  ##########################


lamb_d50 <- seq(1e-3, 0.05, length = 1e2)
pdf(file = "ani_plot_d50.pdf", height = 8, width = 10)
par(mfrow = c(2,2))

plot(lamb_d50, output[[4]][[2]], type = 'l', xlab = "lambda", ylab = "imp_ess", 
     main = "ESS for different lambda")
lines(abline(h=c(0.4,0.6), col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))
lines(abline(v=c(max(lamb_d50[which(output[[4]][[2]] >= 0.4)]), 
    min(lamb_d50[which(output[[4]][[2]] <= 0.6)])), col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb_d50, output[[4]][[3]], type = 'l', xlab = "lambda", ylab = "asymptotic variance", 
     main = "Asymptotic variance for pi^{lambda}")
lines(abline(v=c(max(lamb_d50[which(output[[4]][[2]] >= 0.4)]), 
   min(lamb_d50[which(output[[4]][[2]] <= 0.6)])),  col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))

plot(lamb_d50, output[[4]][[1]], type = 'l', xlab = "lambda", ylab = "asymptotic variance", 
     main = "Asymptotic variance for pi based on IS estimator")
lines(abline(v=c(max(lamb_d50[which(output[[4]][[2]] >= 0.4)]),
   min(lamb_d50[which(output[[4]][[2]] <= 0.6)])),  col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1)))
dev.off()

