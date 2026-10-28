
rm(list = ls())
source("Rcpp/Poisson_functions.R")

##################  Boxplots of marginal efficiency ##################

load("output_poisson.Rdata")

dim <- 51

#### Marginal variance comparison MALA
#x <- c(1:100)
margvar_ism <- matrix(0, nrow = 100, ncol = dim)
margvar_pxm <- matrix(0, nrow = 100, ncol = dim)
for (i in 1:100)
{
  asympmat_ism <- matrix(unlist(output_poisson[[i]][[1]]), nrow = dim, ncol = dim, byrow = T)
  asympmat_pxm <- matrix(unlist(output_poisson[[i]][[2]]), nrow = dim, ncol = dim, byrow = T)
  for (j in 1:dim)
  {
    margvar_ism[i,j] <- asympmat_ism[j,j]
    margvar_pxm[i,j] <- asympmat_pxm[j,j]
  }
}
time_ism <- sapply(output_poisson, function(x) x[[8]])
time_pxm <- sapply(output_poisson, function(x) x[[9]])  

rel_eff_mat_mala <- (margvar_pxm*time_pxm)/(margvar_ism*time_ism)

#### Marginal variance comparison for Barker

margvar_trub <- matrix(0, nrow = 100, ncol = dim)
for (i in 1:100)
{
  asympmat_trub <- matrix(unlist(output_poisson[[i]][[5]]), nrow = dim, ncol = dim, byrow = T)
  for (j in 1:dim)
  {
    margvar_trub[i,j] <- asympmat_trub[j,j]
  }
}

time_trub <- sapply(output_poisson, function(x) x[[12]])
rel_eff_mat_trubark <- (margvar_trub*time_trub)/(margvar_ism*time_ism)

#### Marginal variance comparison for HMC

margvar_ish <- matrix(0, nrow = 100, ncol = dim)
margvar_pxh <- matrix(0, nrow = 100, ncol = dim)
for (i in 1:100)
{
  asympmat_ish <- matrix(unlist(output_poisson[[i]][6]), nrow = dim, ncol = dim, byrow = T)
  asympmat_pxh <- matrix(unlist(output_poisson[[i]][7]), nrow = dim, ncol = dim, byrow = T)
  for (j in 1:dim)
  {
    margvar_ish[i,j] <- asympmat_ish[j,j]
    margvar_pxh[i,j] <- asympmat_pxh[j,j]
  }
}
time_ish <- sapply(output_poisson, function(x) x[[13]])
time_pxh <- sapply(output_poisson, function(x) x[[14]])  

rel_eff_mat_hmc <- (margvar_pxh*time_pxh)/(margvar_ish*time_ish)


# MALA
lower_mala <- apply(rel_eff_mat_mala, 2, quantile, probs = 0.025, na.rm = TRUE)
upper_mala <- apply(rel_eff_mat_mala, 2, quantile, probs = 0.925, na.rm = TRUE)
mu <- colMeans(rel_eff_mat_mala, na.rm = TRUE)

pdf(file = "Rcpp/plots/boxplots_unc_mala.pdf", width = 10, height = 5)
# Plot the band
matplot(1:ncol(rel_eff_mat_mala),cbind(lower_mala, upper_mala, mu),type = "l",lty = c(2, 2, 1), 
        lwd = 2,col = c("blue", "blue", "red"),xlab = "Column",ylab = "Value",
        main = "90% Quantile Band and Mean for MALA")
dev.off()

# Barker
lower_bark <- apply(rel_eff_mat_trubark, 2, quantile, probs = 0.025, na.rm = TRUE)
upper_bark <- apply(rel_eff_mat_trubark, 2, quantile, probs = 0.925, na.rm = TRUE)
mu <- colMeans(rel_eff_mat_trubark, na.rm = TRUE)

pdf(file = "Rcpp/plots/boxplots_unc_bark.pdf", width = 10, height = 5)
# Plot the band
matplot(1:ncol(rel_eff_mat_trubark),cbind(lower_bark, upper_bark, mu),type = "l",lty = c(2, 2, 1), 
        lwd = 2,col = c("blue", "blue", "red"),xlab = "Column",ylab = "Value",
        main = "90% Quantile Band and Mean for Barker")
dev.off() 

# HMC
lower_hmc <- apply(rel_eff_mat_hmc, 2, quantile, probs = 0.025, na.rm = TRUE)
upper_hmc <- apply(rel_eff_mat_hmc, 2, quantile, probs = 0.925, na.rm = TRUE)
mu <- colMeans(rel_eff_mat_hmc, na.rm = TRUE)

pdf(file = "Rcpp/plots/boxplots_unc_hmc.pdf", width = 10, height = 5)
# Plot the band
matplot(1:ncol(rel_eff_mat_hmc),cbind(lower_hmc, upper_hmc, mu),type = "l",lty = c(2, 2, 1), 
        lwd = 2,col = c("blue", "blue", "red"),xlab = "Column",ylab = "Value",
        main = "90% Quantile Band and Mean for HMC")
dev.off() 
