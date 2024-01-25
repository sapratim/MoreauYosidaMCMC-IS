
source("IS_trendf_functions_Pereyra.R")
library("ggplot2")
library("mcmcse")
load("output.Rdata")
lamb_coeff <- 0.001

#  Quantile visualisation
pdf(file = "quantile_plot.pdf")
run_iter <- 20
post_mean <- output[[run_iter]][[1]]
post_med <- output[[run_iter]][[2]]
upper_quant <- output[[run_iter]][[5]]
lower_quant <- output[[run_iter]][[6]]
dataset <- data.frame(x, y, lower_quant, upper_quant, post_med)
plot <- ggplot(dataset, aes(x, y,group = )) + geom_point() +
                        geom_line(aes(x=c(1:100), y=post_med), col = "red")
conf_bands <- plot + geom_ribbon(aes(ymin = lower_quant, ymax = upper_quant), alpha = 0.3) +
                    ggtitle("Piecewise linear model") + labs(x = "argument") + labs(y = "data")
conf_bands
dev.off()
# plot(y, col = "black", main = "95 % credible interval")
# lines(post_mean, type = "l", col = "black")
# lines(post_med, type = "l", col = "yellow")
# 
# lines(upper_quant, type = "l", col = "green")
# lines(lower_quant, type = "l", col = "orange")
# legend("topright", c("observed data", "upper quantile", "lower quantile"), lty = 1,
#        col = c("red", "green", "orange"), cex = 0.8, bty = "n")
# 
