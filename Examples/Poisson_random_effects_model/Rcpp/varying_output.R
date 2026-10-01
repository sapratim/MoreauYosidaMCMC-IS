

############## varying lambda output ################

lambda <- seq(0.0001, 0.004, length = 10)
load("Rcpp/output_variable_lambda.Rdata")

matrices_mala <- lapply(output_lambda, '[[', 1)

mean_mat_mala <- Reduce("+", matrices_mala)/length(matrices_mala)

matrices_hmc <- lapply(output_lambda, '[[', 3)

mean_mat_hmc <- Reduce("+", matrices_hmc)/length(matrices_hmc)

matrices_bark <- lapply(output_lambda, '[[', 2)

mean_mat_bark <- Reduce("+", matrices_bark)/length(matrices_bark)

y_mala <- rowMeans(mean_mat_mala)
y_bark <- rowMeans(mean_mat_bark)
y_hmc  <- rowMeans(mean_mat_hmc)

pdf("Rcpp/plots/varying_lambda.pdf", width = 8, height = 6)

plot(lambda, y_mala,
     type = "o",
     col = "red",
     xlab = expression(lambda),
     ylab = "Relative efficiency per unit time",
     ylim = range(c(y_mala, y_bark, y_hmc)),
     cex.axis = 1.5, cex.lab = 1.5)

lines(lambda, y_bark,
      col = "blue", type = "o")

lines(lambda, y_hmc,
      col = "orange", type = "o")

abline(v = 0.001, lty = 2)

legend("topright",
       legend = c("ISMALA vs PxMALA",
                  "ISMALA vs Barker",
                  "ISHMC vs PxHMC"),
       col = c("red", "blue", "orange"),
       lty = 1)

dev.off()