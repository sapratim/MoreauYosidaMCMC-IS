##########################################
# File makes the plots of interest to us
##########################################
load("Timings_check/small_run.Rdata")
# above loads my and px
# each sampler has samples, acceptance, ADMM counts, final proposals,
# and ADMM counts/roughness at every leapfrog proximal evaluation.


roughness <- function(beta) sum(abs(diff(beta, differences = 2)))

# Existing final-proposal summaries are retained.
my_proposals <- data.frame(iteration = 1:iter, accepted = my$accepted,
                           admm_iterations = my$admm_iterations,
                           roughness = apply(my$proposals, 1, roughness))
px_proposals <- data.frame(iteration = 1:iter, accepted = px$accepted,
                           admm_iterations = px$admm_iterations,
                           roughness = apply(px$proposals, 1, roughness))

# All proximal evaluations (including intermediate leapfrog positions).
# Rows are MCMC iterations, columns are leapfrog steps; unused steps are NA.
leapfrog_data <- function(fit) {
  idx <- which(!is.na(fit$leapfrog_admm), arr.ind = TRUE)
  data.frame(iteration = idx[, 1], leapfrog_step = idx[, 2],
             accepted = fit$accepted[idx[, 1]],
             admm_iterations = fit$leapfrog_admm[idx],
             roughness = fit$leapfrog_roughness[idx])
}
my_leapfrog <- leapfrog_data(my)
px_leapfrog <- leapfrog_data(px)

# calculate roughness of the samples themselves, which may include repeats of the same state
chain_roughness <- data.frame(iteration = 1:iter,
                              MY = apply(my$samples, 1, roughness), PX = apply(px$samples, 1, roughness))

# Roughness versus ADMM iterations, now using ALL leapfrog proximal evaluations.
pdf("roughness_vs_iterations.pdf", height = 5, width = 5)
plot(my_proposals$admm_iterations, my_proposals$roughness, pch = 16, cex = .8, 
     col = adjustcolor("steelblue", alpha.f = .3), xlab = "ADMM Iterations", ylab = "Roughness")
points(px_proposals$admm_iterations, px_proposals$roughness, 
       col = adjustcolor("tomato", alpha.f = .3), pch = 16, cex = .8)
abline(lm(my_proposals$roughness~my_proposals$admm_iterations), lwd = 2, col = "steelblue")
abline(lm(px_proposals$roughness~px_proposals$admm_iterations), lwd = 2, col = "tomato")
legend("bottomleft", legend = c("MY", "PX"), bty = "n",col = c("steelblue", "tomato"), pch = 16)
dev.off()

pdf("admm_iterations.pdf", height = 5, width = 5)
hist(my_leapfrog$admm_iterations, breaks = 20, main = "ADMM Iterations",
     xlab = "ADMM iterations at leapfrog evaluations", ylab = "Counts",
     col = adjustcolor("steelblue", alpha.f = 0.5))
hist(px_leapfrog$admm_iterations, breaks = 20, add = TRUE,
     col = adjustcolor("tomato", alpha.f = 0.5))
legend("topleft", legend = c("MY", "PX"), col = c("steelblue", "tomato"), pch = 15, pt.cex = 2, bty = "n")
dev.off()

pdf("roughness.pdf", height = 5, width = 5)
hist(my_leapfrog$roughness, breaks = 30, main = "Roughness",
     xlab = "Sum of absolute second differences", ylab = "Counts",
     col = adjustcolor("steelblue", alpha.f = 0.5), ylim = c(0, 2.4e4))
hist(px_leapfrog$roughness, breaks = 20, add = TRUE,
     col = adjustcolor("tomato", alpha.f = 0.5))
legend("topleft", legend = c("MY", "PX"), col = c("steelblue", "tomato"), pch = 15, pt.cex = 2, bty = "n")
dev.off()
