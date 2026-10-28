
################################################################################
################## Trendfiltering example output visualisation #################
################################################################################
rm(list = ls())
library(ggplot2)
source("TF_functions.R")

##################  Boxplots of marginal efficiency ##################
load("output_mala.Rdata")
# load("output_bark.Rdata")
load("output_hmc.Rdata")

#### Marginal variance comparison for IS vs PxMALA
x <- c(1:100)
dim <- 100 
margvar_ism <- matrix(0, nrow = 100, ncol = dim)
margvar_pxm <- matrix(0, nrow = 100, ncol = dim)
for (i in 1:100)
{
  asympmat_ism <- matrix(unlist(output_mala[[i]][[3]]), nrow = dim, ncol = dim, byrow = T)
  asympmat_pxm <- matrix(unlist(output_mala[[i]][[4]]), nrow = dim, ncol = dim, byrow = T)
  for (j in 1:dim)
  {
    margvar_ism[i,j] <- asympmat_ism[j,j]
    margvar_pxm[i,j] <- asympmat_pxm[j,j]
  }
}
time_ism <- sapply(output_mala, function(x) x[[5]])
time_pxm <- sapply(output_mala, function(x) x[[6]])
cbind(mean(time_ism), mean(time_pxm))

rel_eff_mat_mala <- (margvar_pxm*time_pxm)/(margvar_ism*time_ism)

#### Marginal variance comparison for IS vs PxHMC
dim <- 100 
margvar_ish <- matrix(0, nrow = 100, ncol = dim)
margvar_pxh <- matrix(0, nrow = 100, ncol = dim)
for (i in 1:100)
{
  asympmat_ish <- matrix(unlist(output_hmc[[i]][3]), nrow = dim, ncol = dim, byrow = T)
  asympmat_pxh <- matrix(unlist(output_hmc[[i]][4]), nrow = dim, ncol = dim, byrow = T)
  for (j in 1:dim)
  {
    margvar_ish[i,j] <- asympmat_ish[j,j]
    margvar_pxh[i,j] <- asympmat_pxh[j,j]
  }
}

time_ish <- sapply(output_hmc, function(x) x[[5]])
time_pxh <- sapply(output_hmc, function(x) x[[6]])
cbind(mean(time_ish), mean(time_pxh))

rel_eff_mat_hmc <- (margvar_pxh*time_pxh)/(margvar_ish*time_ish)

###################. Boxplots ####################

boxplot(rel_eff_mat_mala, main = "Column-wise Boxplots", xlab = "Columns",
    ylab = "Values", col = "lightblue",names = paste0("Column ", 1:ncol(rel_eff_mat_mala)))
boxplot(rel_eff_mat_hmc, main = "Column-wise Boxplots", xlab = "Columns",
        ylab = "Values", col = "lightblue",names = paste0("Column ", 1:ncol(rel_eff_mat_mala)))


# MALA
lower_mala <- apply(rel_eff_mat_mala, 2, quantile, probs = 0.025, na.rm = TRUE)
upper_mala <- apply(rel_eff_mat_mala, 2, quantile, probs = 0.925, na.rm = TRUE)
mu <- colMeans(rel_eff_mat_mala, na.rm = TRUE)

pdf(file = "plots/boxplots_unc_mala.pdf", width = 10, height = 5)
# Plot the band
matplot(1:ncol(rel_eff_mat_mala),cbind(lower_mala, upper_mala, mu),type = "l",lty = c(2, 2, 1), 
                  lwd = 2,col = c("blue", "blue", "red"),xlab = "Column",ylab = "Value",
                                                      main = "90% Quantile Band and Mean for MALA")
dev.off()

# HMC
lower_hmc <- apply(rel_eff_mat_hmc, 2, quantile, probs = 0.025, na.rm = TRUE)
upper_hmc <- apply(rel_eff_mat_hmc, 2, quantile, probs = 0.925, na.rm = TRUE)
mu <- colMeans(rel_eff_mat_hmc, na.rm = TRUE)

pdf(file = "plots/boxplots_unc_hmc.pdf", width = 10, height = 5)
# Plot the band
matplot(1:ncol(rel_eff_mat_hmc),cbind(lower_hmc, upper_hmc, mu),type = "l",lty = c(2, 2, 1), 
        lwd = 2,col = c("blue", "blue", "red"),xlab = "Column",ylab = "Value",
        main = "90% Quantile Band and Mean for HMC")
dev.off() 





