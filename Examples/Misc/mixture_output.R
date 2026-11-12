###########################################
# mixture_output.R
# Loads the list saved by mixture_run.R and makes the plots.
# No simulation is run here.
###########################################

out_dir <- "mixture_results"

results <- readRDS(file.path(out_dir, "mixture_results.rds"))

sigma2_vec <- results$settings$sigma2_vec

# 1. Laplace / Gaussian

laplace_gaussian_mean_ess <- results$laplace_gaussian$mean_ess
laplace_gaussian_sigma2 <- results$laplace_gaussian$sigma2
laplace_gaussian_mix_logs <- results$laplace_gaussian$mix_logs
laplace_gaussian_my_logs <- results$laplace_gaussian$my_logs

pdf(file.path(out_dir, "laplace_gaussian_ess.pdf"), width = 10, height = 6)

par(mar = c(5.5, 5.5, 4, 2) + 0.1)
plot(
  sigma2_vec, laplace_gaussian_mean_ess,
  type = "b", pch = 19, col = "steelblue",
  xlab = expression(sigma^2), ylab = expression("Average " * n[e] / n),
  cex.lab = 2, cex.axis = 2
)
abline(v = laplace_gaussian_sigma2, col = "firebrick", lty = 2)

dev.off()

laplace_gaussian_limits <- range(
  laplace_gaussian_mix_logs, laplace_gaussian_my_logs
)
laplace_gaussian_breaks <- pretty(laplace_gaussian_limits, n = 300)

laplace_gaussian_mix_hist <- hist(
  laplace_gaussian_mix_logs,
  breaks = laplace_gaussian_breaks, plot = FALSE
)
laplace_gaussian_my_hist <- hist(
  laplace_gaussian_my_logs,
  breaks = laplace_gaussian_breaks, plot = FALSE
)

laplace_gaussian_ylim <- c(
  0, max(laplace_gaussian_mix_hist$density, laplace_gaussian_my_hist$density)
)

pdf(file.path(out_dir, "laplace_gaussian_weights.pdf"), width = 10, height = 6)

par(mar = c(5.5, 5.5, 4, 2) + 0.1)
hist(
  laplace_gaussian_my_logs, breaks = laplace_gaussian_breaks, probability = TRUE,
  xlim = laplace_gaussian_limits, ylim = laplace_gaussian_ylim,
  col = rgb(1, 0, 0, 0.4), border = "red",
  main = "Histogram of Laplace",
  xlab = "Log importance weight", ylab = "Density",
  cex.lab = 2, cex.axis = 2, cex.main = 1.8
)
hist(
  laplace_gaussian_mix_logs, breaks = laplace_gaussian_breaks, probability = TRUE,
  col = rgb(0, 0, 1, 0.4), border = "blue", add = TRUE
)
legend(
  "topright", c("MY-MALA IS", "Mixture IS"), bty = "n", cex = 1.5,
  fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4))
)

dev.off()

laplace_gaussian_summary <- results$laplace_gaussian$summary

# 2. Laplace / Student t

laplace_student_mean_ess <- results$laplace_student$mean_ess
laplace_student_sigma2 <- results$laplace_student$sigma2
laplace_student_mix_logs <- results$laplace_student$mix_logs
laplace_student_my_logs <- results$laplace_student$my_logs

pdf(file.path(out_dir, "laplace_student_ess.pdf"), width = 10, height = 6)

par(mar = c(5.5, 5.5, 4, 2) + 0.1)
plot(
  sigma2_vec, laplace_student_mean_ess,
  type = "b", pch = 19, col = "steelblue",
  xlab = expression(sigma^2), ylab = expression("Average " * n[e] / n),
  cex.lab = 2, cex.axis = 2
)
abline(v = laplace_student_sigma2, col = "firebrick", lty = 2)

dev.off()

laplace_student_limits <- range(
  laplace_student_mix_logs, laplace_student_my_logs
)
laplace_student_breaks <- pretty(laplace_student_limits, n = 300)

laplace_student_mix_hist <- hist(
  laplace_student_mix_logs,
  breaks = laplace_student_breaks, plot = FALSE
)
laplace_student_my_hist <- hist(
  laplace_student_my_logs,
  breaks = laplace_student_breaks, plot = FALSE
)

laplace_student_ylim <- c(
  0, max(laplace_student_mix_hist$density, laplace_student_my_hist$density)
)

pdf(file.path(out_dir, "laplace_student_weights.pdf"), width = 10, height = 6)

par(mar = c(5.5, 5.5, 4, 2) + 0.1)
hist(
  laplace_student_my_logs, breaks = laplace_student_breaks, probability = TRUE,
  xlim = laplace_student_limits, ylim = laplace_student_ylim,
  col = rgb(1, 0, 0, 0.4), border = "red",
  main = "Histogram of Laplace",
  xlab = "Log importance weight", ylab = "Density",
  cex.lab = 2, cex.axis = 2, cex.main = 1.8
)
hist(
  laplace_student_mix_logs, breaks = laplace_student_breaks, probability = TRUE,
  col = rgb(0, 0, 1, 0.4), border = "blue", add = TRUE
)
legend(
  "topright", c("MY-MALA IS", "Mixture IS"), bty = "n", cex = 1.5,
  fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4))
)

dev.off()

laplace_student_summary <- results$laplace_student$summary

# 3. Quartic / Gaussian

exp4_gaussian_mean_ess <- results$exp4_gaussian$mean_ess
exp4_gaussian_sigma2 <- results$exp4_gaussian$sigma2
exp4_gaussian_mix_logs <- results$exp4_gaussian$mix_logs
exp4_gaussian_my_logs <- results$exp4_gaussian$my_logs

pdf(file.path(out_dir, "exp4_gaussian_ess.pdf"), width = 10, height = 6)

par(mar = c(5.5, 5.5, 4, 2) + 0.1)
plot(
  sigma2_vec, exp4_gaussian_mean_ess,
  type = "b", pch = 19, col = "steelblue",
  xlab = expression(sigma^2), ylab = expression("Average " * n[e] / n),
  cex.lab = 2, cex.axis = 2
)
abline(v = exp4_gaussian_sigma2, col = "firebrick", lty = 2)

dev.off()

exp4_gaussian_limits <- range(
  exp4_gaussian_mix_logs, exp4_gaussian_my_logs
)
exp4_gaussian_breaks <- pretty(exp4_gaussian_limits, n = 300)

exp4_gaussian_mix_hist <- hist(
  exp4_gaussian_mix_logs,
  breaks = exp4_gaussian_breaks, plot = FALSE
)
exp4_gaussian_my_hist <- hist(
  exp4_gaussian_my_logs,
  breaks = exp4_gaussian_breaks, plot = FALSE
)

exp4_gaussian_ylim <- c(
  0, max(exp4_gaussian_mix_hist$density, exp4_gaussian_my_hist$density)
)

pdf(file.path(out_dir, "exp4_gaussian_weights.pdf"), width = 10, height = 6)

par(mar = c(5.5, 5.5, 4, 2) + 0.1)
hist(
  exp4_gaussian_my_logs, breaks = exp4_gaussian_breaks, probability = TRUE,
  xlim = exp4_gaussian_limits, ylim = exp4_gaussian_ylim,
  col = rgb(1, 0, 0, 0.4), border = "red",
  main = "Histogram of Super-Gaussian",
  xlab = "Log importance weight", ylab = "Density",
  cex.lab = 2, cex.axis = 2, cex.main = 1.8
)
hist(
  exp4_gaussian_mix_logs, breaks = exp4_gaussian_breaks, probability = TRUE,
  col = rgb(0, 0, 1, 0.4), border = "blue", add = TRUE
)
legend(
  "topleft", c("MY-MALA IS", "Mixture IS"), bty = "n", cex = 1.5,
  fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4))
)

dev.off()

exp4_gaussian_summary <- results$exp4_gaussian$summary

# 4. Quartic / Student t

exp4_student_mean_ess <- results$exp4_student$mean_ess
exp4_student_sigma2 <- results$exp4_student$sigma2
exp4_student_mix_logs <- results$exp4_student$mix_logs
exp4_student_my_logs <- results$exp4_student$my_logs

pdf(file.path(out_dir, "exp4_student_ess.pdf"), width = 10, height = 6)

par(mar = c(5.5, 5.5, 4, 2) + 0.1)
plot(
  sigma2_vec, exp4_student_mean_ess,
  type = "b", pch = 19, col = "steelblue",
  xlab = expression(sigma^2), ylab = expression("Average " * n[e] / n),
  cex.lab = 2, cex.axis = 2
)
abline(v = exp4_student_sigma2, col = "firebrick", lty = 2)

dev.off()

exp4_student_limits <- range(
  exp4_student_mix_logs, exp4_student_my_logs
)
exp4_student_breaks <- pretty(exp4_student_limits, n = 300)

exp4_student_mix_hist <- hist(
  exp4_student_mix_logs,
  breaks = exp4_student_breaks, plot = FALSE
)
exp4_student_my_hist <- hist(
  exp4_student_my_logs,
  breaks = exp4_student_breaks, plot = FALSE
)

exp4_student_ylim <- c(
  0, max(exp4_student_mix_hist$density, exp4_student_my_hist$density)
)

pdf(file.path(out_dir, "exp4_student_weights.pdf"), width = 10, height = 6)

par(mar = c(5.5, 5.5, 4, 2) + 0.1)
hist(
  exp4_student_my_logs, breaks = exp4_student_breaks, probability = TRUE,
  xlim = exp4_student_limits, ylim = exp4_student_ylim,
  col = rgb(1, 0, 0, 0.4), border = "red",
  main = "Histogram of Super-Gaussian",
  xlab = "Log importance weight", ylab = "Density",
  cex.lab = 2, cex.axis = 2, cex.main = 1.8
)
hist(
  exp4_student_mix_logs, breaks = exp4_student_breaks, probability = TRUE,
  col = rgb(0, 0, 1, 0.4), border = "blue", add = TRUE
)
legend(
  "topleft", c("MY-MALA IS", "Mixture IS"), bty = "n", cex = 1.5,
  fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4))
)

dev.off()

exp4_student_summary <- results$exp4_student$summary
