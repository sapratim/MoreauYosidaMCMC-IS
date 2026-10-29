
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







