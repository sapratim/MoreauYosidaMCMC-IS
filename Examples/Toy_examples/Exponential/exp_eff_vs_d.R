
#####################################################################
################# Exponential output visualisation ##################
#####################################################################
set.seed(100)
library(mcmcse)
source("exp_d4_functions.R")
load("output_d1.Rdata")
load("output_d5.Rdata")
load("output_d10.Rdata")
load("output_d20.Rdata")
load("output_d50.Rdata")
load("output_d100.Rdata")

eff_mat_d1 <- sapply(output_exp_d1, function(l) l[[3]])
eff_mat_d5 <- colMeans(sapply(output_exp_d5, function(l) l[[3]]))
eff_mat_d10 <- colMeans(sapply(output_exp_d10, function(l) l[[3]]))
eff_mat_d20 <- colMeans(sapply(output_exp_d20, function(l) l[[3]]))
eff_mat_d50 <- colMeans(sapply(output_exp_d50, function(l) l[[3]]))
eff_mat_d100 <- colMeans(sapply(output_exp_d100, function(l) l[[3]]))

se_d1 <- sd(eff_mat_d1)/sqrt(length(eff_mat_d1))
se_d5 <- sd(eff_mat_d5)/sqrt(length(eff_mat_d5))
se_d10 <- sd(eff_mat_d10)/sqrt(length(eff_mat_d10))
se_d20 <- sd(eff_mat_d20)/sqrt(length(eff_mat_d20))
se_d50 <- sd(eff_mat_d50)/sqrt(length(eff_mat_d50))
se_d100 <- sd(eff_mat_d100)/sqrt(length(eff_mat_d100))

######## d vs Relative efficiency ########

mean_vec <- c(mean(eff_mat_d1), mean(eff_mat_d5), mean(eff_mat_d10), mean(eff_mat_d20),
              mean(eff_mat_d50), mean(eff_mat_d100))
se_vec <- c(se_d1, se_d5, se_d10, se_d20, se_d50, se_d100)

dim <- c(1,5,10,20,50,100)

pdf(file = "plots/dimension_vs_re.pdf", width = 8, height = 6)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(dim, mean_vec, type = 'o', main = "Super Gaussian",
     ylim = c(0.5, 2), xlab = "dimension", ylab = "Relative efficiency",
     cex.lab = 1.8,  cex.axis = 1.5,  cex.main = 1.8)
segments(dim, mean_vec - se_vec,
         dim, mean_vec + se_vec, lwd = 3, col = "blue")
dev.off()