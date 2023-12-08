
load("mala_chain.Rdata")
load("weights.Rdata")

#   Difference in modes
freq_mode <- trendfilter(x,y, k=1,lambda = sigma2_hat*alpha_hat,
                         control = trendfilter.control.list(obj_tol = tol, max_iter = 1e4L))$beta
proxval <- prox_func(freq_mode, lamb_coeff,alpha_hat, sigma2_hat, k=1, grid = x)
mode_diff <- abs(freq_mode-proxval)
mode_diff   

i <- 1

# Repeated execution gives density functions for different components
plot(density(mala.is[[1]][, i]))
abline(v = proxval[i], col = "red")
abline(v= freq_mode[i], col = "blue")
legend("topright", c("est_density_MCMC", "frequentist mode", "prox value at mode"), lty = 1,
       col = c("black", "blue", "red"), cex = 0.6, bty = "n")
i <- i + 1

###  Trace plots
plot.ts(mala.is[[1]][, 1:10])
plot.ts(mala.is[[1]][, 11:20])
plot.ts(mala.is[[1]][, 21:30])
plot.ts(mala.is[[1]][, 31:40])
plot.ts(mala.is[[1]][, 41:50])
plot.ts(mala.is[[1]][, 51:60])
plot.ts(mala.is[[1]][, 61:70])
plot.ts(mala.is[[1]][, 71:80])
plot.ts(mala.is[[1]][, 81:90])
plot.ts(mala.is[[1]][, 91:100])