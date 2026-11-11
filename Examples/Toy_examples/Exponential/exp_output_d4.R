# Run exp_run_d4.R first: these plots use replication means.
# *_reps matrices have lambda values in rows and replications in columns.

# Pointwise 95% t confidence intervals for the mean across replications.
plot_mean_ci <- function(x, rep_values, ...) {
        n <- ncol(rep_values)
        stopifnot(n >= 2)
        avg <- rowMeans(rep_values)
        margin <- qt(0.975, df = n - 1) * apply(rep_values, 1, sd) / sqrt(n)
        lower <- avg - margin
        upper <- avg + margin
        plot(x, avg, type = "n", ylim = range(lower, upper), ...)
        polygon(c(x, rev(x)), c(lower, rev(upper)),
                col = adjustcolor("grey50", alpha.f = 0.25), border = NA)
        lines(x, avg)
}

source("exp_d4_functions.R")
load("output_d1.Rdata")
load("output_d5.Rdata")
load("output_d10.Rdata")
load("output_d20.Rdata")

iter <- 1e4    ##### Number of iterations
n_lambda <- length(output_exp_d1[[1]])
reps <- length(output_exp_d1)
cutoff <- c(0.4,0.8)
lamb_d1 <- seq(1e-4, 5, length = n_lambda)
lamb_d5 <- seq(1e-4, .6, length = n_lambda)
lamb_d10 <- seq(1e-3, 0.5, length = n_lambda)
lamb_d20 <- seq(1e-3, 0.2, length = n_lambda)


##########################  For d = 1  ##########################
load("output_d1.Rdata")

ess_d1_reps <- vapply(output_exp_d1, function(rep_output)
        vapply(rep_output, function(l) mean(l[[4]]), numeric(1)),
        numeric(length(lamb_d1)))
ess_d1 <- rowMeans(ess_d1_reps)
pilam_var_d1_reps <- vapply(output_exp_d1, function(rep_output)
        vapply(rep_output, function(l) mean(l[[2]]), numeric(1)),
        numeric(length(lamb_d1)))
pilam_var_d1 <- rowMeans(pilam_var_d1_reps)
is_var_d1_reps <- vapply(output_exp_d1, function(rep_output)
        vapply(rep_output, function(l) mean(l[[1]]), numeric(1)),
        numeric(length(lamb_d1)))
is_var_d1 <- rowMeans(is_var_d1_reps)

opt_win <- c(max(lamb_d1[which(ess_d1 >= cutoff[1])]), min(lamb_d1[which(ess_d1 <= cutoff[2])]))


pdf(file = "plots/exp_plot_d1.pdf", height = 4, width = 12)
par(mar = c(5, 7, 4, 2), cex.lab = 1.5, cex.axis = 1.5)
par(mfrow = c(1,3))

plot_mean_ci(lamb_d1, ess_d1_reps, xlab = expression(lambda), ylab = expression(n[e]/n))
# lines(abline(h=cutoff, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))
lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))

plot_mean_ci(lamb_d1, pilam_var_d1_reps, xlab = expression(lambda), 
             ylab = expression(paste("MCMC ESS under ",pi^lambda))) 
# lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))

plot_mean_ci(lamb_d1, is_var_d1_reps, xlab = expression(lambda),
             ylab ="Importance sampling asymptotic variance")
lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))
dev.off()

##########################  For d = 5  ##########################
load("output_d5.Rdata")

ess_d5_reps <- vapply(output_exp_d5, function(rep_output)
        vapply(rep_output, function(l) mean(l[[4]]), numeric(1)),
        numeric(length(lamb_d5)))
ess_d5 <- rowMeans(ess_d5_reps)
pilam_var_d5_reps <- vapply(output_exp_d5, function(rep_output)
        vapply(rep_output, function(l) mean(l[[2]]), numeric(1)),
        numeric(length(lamb_d5)))
pilam_var_d5 <- rowMeans(pilam_var_d5_reps)
is_var_d5_reps <- vapply(output_exp_d5, function(rep_output)
        vapply(rep_output, function(l) mean(l[[1]]), numeric(1)),
        numeric(length(lamb_d5)))
is_var_d5 <- rowMeans(is_var_d5_reps)

opt_win <- c(max(lamb_d5[which(ess_d5 >= cutoff[1])]), min(lamb_d5[which(ess_d5 <= cutoff[2])]))


pdf(file = "plots/exp_plot_d5.pdf", height = 4, width = 12)
par(mar = c(5, 7, 4, 2), cex.lab = 1.5, cex.axis = 1.5)
par(mfrow = c(1,3))

plot_mean_ci(lamb_d5, ess_d5_reps, xlab = expression(lambda), ylab = expression(n[e]/n))
# lines(abline(h=cutoff, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))
lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))

plot_mean_ci(lamb_d5, pilam_var_d5_reps, xlab = expression(lambda), 
             ylab = expression(paste("MCMC ESS under ",pi^lambda))) 
# lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))

plot_mean_ci(lamb_d5, is_var_d5_reps, xlab = expression(lambda), 
             ylab ="Importance sampling asymptotic variance")

lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))
dev.off()


##########################  For d = 10  ##########################
load("output_d10.Rdata")

ess_d10_reps <- vapply(output_exp_d10, function(rep_output)
        vapply(rep_output, function(l) mean(l[[4]]), numeric(1)),
        numeric(length(lamb_d10)))
ess_d10 <- rowMeans(ess_d10_reps)
pilam_var_d10_reps <- vapply(output_exp_d10, function(rep_output)
        vapply(rep_output, function(l) mean(l[[2]]), numeric(1)),
        numeric(length(lamb_d10)))
pilam_var_d10 <- rowMeans(pilam_var_d10_reps)
is_var_d10_reps <- vapply(output_exp_d10, function(rep_output)
        vapply(rep_output, function(l) mean(l[[1]]), numeric(1)),
        numeric(length(lamb_d10)))
is_var_d10 <- rowMeans(is_var_d10_reps)

opt_win <- c(max(lamb_d10[which(ess_d10 >= cutoff[1])]), min(lamb_d10[which(ess_d10 <= cutoff[2])]))

pdf(file = "plots/exp_plot_d10.pdf", height = 4, width = 12)
par(mar = c(5, 7, 4, 2), cex.lab = 1.5, cex.axis = 1.5)
par(mfrow = c(1,3))

plot_mean_ci(lamb_d10, ess_d10_reps, xlab = expression(lambda), ylab = expression(n[e]/n))
# lines(abline(h=cutoff, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))
lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))

plot_mean_ci(lamb_d10, pilam_var_d10_reps, xlab = expression(lambda), 
             ylab = expression(paste("MCMC ESS under ",pi^lambda))) 
# lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))

plot_mean_ci(lamb_d10, is_var_d10_reps, xlab = expression(lambda), 
             ylab ="Importance sampling asymptotic variance")
lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))
dev.off()


##########################  For d = 20  ##########################
load("output_d20.Rdata")

ess_d20_reps <- vapply(output_exp_d20, function(rep_output)
        vapply(rep_output, function(l) mean(l[[4]]), numeric(1)),
        numeric(length(lamb_d20)))
ess_d20 <- rowMeans(ess_d20_reps)
pilam_var_d20_reps <- vapply(output_exp_d20, function(rep_output)
        vapply(rep_output, function(l) mean(l[[2]]), numeric(1)),
        numeric(length(lamb_d20)))
pilam_var_d20 <- rowMeans(pilam_var_d20_reps)
is_var_d20_reps <- vapply(output_exp_d20, function(rep_output)
        vapply(rep_output, function(l) mean(l[[1]]), numeric(1)),
        numeric(length(lamb_d20)))
is_var_d20 <- rowMeans(is_var_d20_reps)

opt_win <- c(max(lamb_d20[which(ess_d20 >= cutoff[1])]), min(lamb_d20[which(ess_d20 <= cutoff[2])]))


pdf(file = "plots/exp_plot_d20.pdf", height = 4, width = 12)
par(mar = c(5, 7, 4, 2), cex.lab = 1.5, cex.axis = 1.5)
par(mfrow = c(1,3))

plot_mean_ci(lamb_d20, ess_d20_reps, xlab = expression(lambda), ylab = expression(n[e]/n))
# lines(abline(h=cutoff, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))
lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))

plot_mean_ci(lamb_d20, pilam_var_d20_reps, xlab = expression(lambda), 
             ylab = expression(paste("MCMC ESS under ",pi^lambda))) 
# lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))

plot_mean_ci(lamb_d20, is_var_d20_reps, xlab = expression(lambda), 
             ylab ="Importance sampling asymptotic variance")
lines(abline(v=opt_win, col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1)))
dev.off()
