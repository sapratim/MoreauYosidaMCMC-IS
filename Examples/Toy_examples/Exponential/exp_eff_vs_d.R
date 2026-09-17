
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

######## d vs Relative efficiency ########

rel_eff_d_1 <- mean(sapply(output_exp_d1, function(l) l[[3]]))
rel_eff_d_5 <- mean(colMeans(sapply(output_exp_d5, function(l) l[[3]])))
rel_eff_d_10 <- mean(colMeans(sapply(output_exp_d10, function(l) l[[3]])))
rel_eff_d_20 <- mean(colMeans(sapply(output_exp_d20, function(l) l[[3]])))
rel_eff_d_100 <- mean(colMeans(sapply(output_exp_d50, function(l) l[[3]])))
rel_eff_d_500 <- mean(colMeans(sapply(output_exp_d100, function(l) l[[3]])))

rel_eff_vec <- c(rel_eff_d_1, rel_eff_d_5, rel_eff_d_10, 
                    rel_eff_d_20, rel_eff_d_100, rel_eff_d_500)
dim <- c(1,5,10,20,50,100)

pdf(file = "plots/dimension_vs_re.pdf", width = 8, height = 6)
plot(dim, rel_eff_vec, type = 'o', main = "Super-Gaussian",
     ylim = c(0.5, 2), xlab = "dimension", ylab = "Relative efficiency")
dev.off()