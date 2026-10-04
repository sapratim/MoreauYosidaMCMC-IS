

############## varying lambda output ################

lambda <- seq(0.0001, 0.004, length = 10)
load("Rcpp/output_variable_lambda.Rdata")

################################ MALA ###############################

matrices_mala <- lapply(output_lambda, '[[', 1)
se_mala_mat <- sapply(matrices_mala, function(x) rowMeans(x)) #avg eff over components for all lambda
mean_val_mala <- rowMeans(se_mala_mat)
se_val_mala <- apply(se_mala_mat, 1, sd) / sqrt(ncol(se_mala_mat))

# lowerquant_mala <- apply(se_mala_mat, 1, quantile, probs = 0.025)
# median_mala <- apply(se_mala_mat, 1, quantile, probs = 0.5)
# upperquant_mala <- apply(se_mala_mat, 1, quantile, probs = 0.975)

#mean_mat_mala <- Reduce("+", matrices_mala)/length(matrices_mala)

################################ HMC #################################

matrices_hmc <- lapply(output_lambda, '[[', 3)
se_hmc_mat <- sapply(matrices_hmc, function(x) rowMeans(x)) #avg eff over components for all lambda
mean_val_hmc <- rowMeans(se_hmc_mat)
se_val_hmc <- apply(se_hmc_mat, 1, sd) / sqrt(ncol(se_hmc_mat))

# lowerquant_hmc <- apply(se_hmc_mat, 1, quantile, probs = 0.025)
# median_hmc <- apply(se_hmc_mat, 1, quantile, probs = 0.5)
# upperquant_hmc <- apply(se_hmc_mat, 1, quantile, probs = 0.975)

#mean_mat_hmc <- Reduce("+", matrices_hmc)/length(matrices_hmc)

################################ Barker #################################
matrices_bark <- lapply(output_lambda, '[[', 2)
se_bark_mat <- sapply(matrices_bark, function(x) rowMeans(x)) #avg eff over components for all lambda
mean_val_bark <- rowMeans(se_bark_mat)
se_val_bark <- apply(se_bark_mat, 1, sd) / sqrt(ncol(se_bark_mat))

# lowerquant_bark <- apply(se_bark_mat, 1, quantile, probs = 0.025)
# median_bark <- apply(se_bark_mat, 1, quantile, probs = 0.5)
# upperquant_bark <- apply(se_bark_mat, 1, quantile, probs = 0.975)

#mean_mat_bark <- Reduce("+", matrices_bark)/length(matrices_bark)

# y_mala <- rowMeans(mean_mat_mala)
# y_bark <- rowMeans(mean_mat_bark)
# y_hmc  <- rowMeans(mean_mat_hmc)

pdf("Rcpp/plots/varying_lambda.pdf", width = 8, height = 6)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(lambda, mean_val_mala,
     type = "l",
     col = "red", log = "y",
     xlab = expression(lambda),
     ylab = "Log average relative efficiency ",
     ylim = range(c(mean_val_mala, mean_val_bark, mean_val_hmc)),
     cex.axis = 1.5, cex.lab = 1.5)
segments(lambda, mean_val_mala - se_val_mala,
         lambda, mean_val_mala + se_val_mala,
         col = "red", lwd = 5)

lines(lambda, mean_val_bark,
      col = "blue", type = "l")
# lines(lambda, lowerquant_bark, col = "blue", lty = 4)
# lines(lambda, upperquant_bark, col = "blue", lty = 4)
segments(lambda, mean_val_bark - se_val_bark,
         lambda, mean_val_bark + se_val_bark,
         col = "blue", lwd = 5)

lines(lambda, mean_val_hmc,
      col = "orange", type = "l")
# lines(lambda, lowerquant_hmc, col = "orange", lty = 4)
# lines(lambda, upperquant_hmc, col = "orange", lty = 4)
segments(lambda, mean_val_hmc - se_val_hmc,
         lambda, mean_val_hmc + se_val_hmc,
         col = "orange", lwd = 5)

abline(v = 0.001, lty = 2)

legend("topright",
       legend = c("ISMALA vs PxMALA",
                  "ISMALA vs Barker",
                  "ISHMC vs PxHMC"),
       col = c("red", "blue", "orange"),
       lty = 1, cex = 0.72)

dev.off()