
source("IS_trendf_functions_Pereyra.R")
load("mala_chain.Rdata")
load("weights.Rdata")
lamb_coeff <- 0.0005
iter <- 1e4

#   Difference in modes
freq_mode <- trendfilter(x,y, k=1,lambda = sigma2_hat*alpha_hat,
                         control = trendfilter.control.list(obj_tol = tol, max_iter = 1e4L))$beta
proxval <- prox_func(freq_mode, lamb_coeff,alpha_hat, sigma2_hat, k=1, grid = x)
mode_diff <- abs(freq_mode-proxval)
mode_diff   

i <- 1

# Repeated execution gives density functions for different components
plot(density(mala_chain[, i]))
abline(v = proxval[i], col = "red")
abline(v= freq_mode[i], col = "blue")
legend("topright", c("est_density_MCMC", "frequentist mode", "prox value at mode"), lty = 1,
       col = c("black", "blue", "red"), cex = 0.6, bty = "n")
i <- i + 1

###  Trace plots
plot.ts(mala_chain[, 1:10])
plot.ts(mala_chain[, 11:20])
plot.ts(mala_chain[, 21:30])
plot.ts(mala_chain[, 31:40])
plot.ts(mala_chain[, 41:50])
plot.ts(mala_chain[, 51:60])
plot.ts(mala_chain[, 61:70])
plot.ts(mala_chain[, 71:80])
plot.ts(mala_chain[, 81:90])
plot.ts(mala_chain[, 91:100])

hist(weights)

weight_mat <- matrix(0, nrow = iter, ncol = length(y))

for (i in 1:iter) {
  weight_mat[i,] <- mala_chain[i,]*exp(weights[i])
}

num_sum <- apply(weight_mat, 2, sum)
weights_sum <- sum(exp(weights))
post_mean <- num_sum/weights_sum

plot(y, type = "o", col = "red", main = "data vs posterior mean")
lines(post_mean, col = "green")
legend("topright", c("observed data", "posterior IS mean"), lty = 1,
       col = c("red", "green"), cex = 0.8, bty = "n")


#  Quantile visualisation

augm_mat <- cbind(mala_chain,weights)

quant <- function(j)
{
  mat_ordered <- augm_mat[order(augm_mat[,j], decreasing = FALSE), ]
  order_comp <- mat_ordered[,j]
  weights_order <- mat_ordered[, length(y)+1]
  bound_mat <- cbind(order_comp, weights_order)
  return(bound_mat)
}

upper_quant <- numeric(length = length(y))
lower_quant <- numeric(length = length(y))
signif_level <- 0.025

for (i in 1:length(y)) 
 {
  initial_mat <- quant(i)
  mat_sum <- apply(initial_mat, 2, sum)
  wts_prop <- initial_mat[,2]/mat_sum[2]
  final_mat <- cbind(initial_mat[1,], cumsum(wts_prop))
  lower_index <- min(which(final_mat[,2] >= signif_level))
  upper_index <- min(which(final_mat[,2] >= (1 - signif_level)))
  upper_quant[i] <- initial_mat[upper_index,1]
  lower_quant[i] <- initial_mat[lower_index,1]
}

dataset <- data.frame(x, y, lower_quant, upper_quant)
library("ggplot2")
plot <- ggplot(dataset, aes(x, y)) + geom_point()
conf_bands <- plot + geom_ribbon(aes(ymin = lower_quant, ymax = upper_quant),  
              alpha = 0.2, fill="orange", color="orange")
conf_bands

# plot(y, type = "l", col = "red", main = "95 % credible interval")
# lines(upper_quant, type = "l", col = "green")
# lines(lower_quant, type = "l", col = "orange")
# legend("topright", c("observed data", "upper quantile", "lower quantile"), lty = 1,
#        col = c("red", "green", "orange"), cex = 0.8, bty = "n")
