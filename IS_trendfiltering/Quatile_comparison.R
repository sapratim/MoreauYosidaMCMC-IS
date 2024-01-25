upper_quant_pxm <- numeric(length = length(y))
lower_quant_pxm <- numeric(length = length(y))
for(i in 1:length(y))
{
  upper_quant_pxm[i] <- quantile(pxmala.run[,i], probs = 0.975)
  lower_quant_pxm[i] <- quantile(pxmala.run[,i], probs = 0.025)
}

pdf(file = "quantile_sample_plot.pdf", width = 7, height = 5)
plot(y, col = "black", main = "95 % credible interval")
lines(post_mean, type = "l", col = "red")
# lines(post_med, type = "l", col = "yellow")
 
lines(upper_quant, type = "o", col = "green")
lines(lower_quant, type = "o", col = "orange")
lines(upper_quant_pxm, type = "o", col = "blue")
lines(lower_quant_pxm, type = "o", col = "purple")
legend("topright", c("posterior mean", "upper quantile_is", "lower quantile_is",
                     "upper quantile_pxm", "lower_quantile_pxm"), lty = 1,
     col = c("red", "green", "orange", "blue", "purple"), cex = 0.6, bty = "n")
dev.off()