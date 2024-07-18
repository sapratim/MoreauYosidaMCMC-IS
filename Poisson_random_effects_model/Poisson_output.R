################################################################################
################## Poisson example output visualisation ########################
################################################################################

rm(list = ls())
source("Poisson_functions.R")
load("output_single.Rdata")
load("output_poisson.Rdata")

################ ACF plots ################
dim <- 51
pdf("plots/acf_poisson.pdf", height = 6, width = 12)
par(mfrow = c(1,3))

### MALA
lag.max <- 50
acf_ism <- acf(output_chain[[1]][,1], plot = FALSE, lag.max = lag.max)$acf
acf_pxm <- acf(output_chain[[2]][,1], plot = FALSE, lag.max = lag.max)$acf

diff.acf <- matrix(0, ncol = dim, nrow = lag.max + 1)
diff.acf[,1] <- acf_ism - acf_pxm

# plot(1:length(acf_ism), acf_ism - acf_pxm, col = "blue", type = "l",
#      xlab = "Lag", ylab = "Autocorrelation", ylim = c(-0.8, .3))
# # lines(1:length(acf_pxm), acf_pxm, col = "red", type = "l")
# legend("bottomright", c("MYMALA", "PxMALA"), lty = 1,
#        col = c("blue", "red"), cex = 0.75, bty = "n")

for (i in 2:dim) 
{
  acf_ism <- acf(output_chain[[1]][,i], plot = FALSE,lag.max = lag.max)$acf
  acf_pxm <- acf(output_chain[[2]][,i], plot = FALSE,lag.max = lag.max)$acf
  diff.acf[,i] <- acf_ism - acf_pxm
  #lines(1:length(acf_ism), acf_ism - acf_pxm, col = "blue", type = "l")
  #lines(1:length(acf_pxm), acf_pxm, col = "red", type = "l")
}

# Make boxplot of ACFs
boxplot(t(diff.acf),
  xlab = "Lags", ylab = "Difference in ACF (MYMala - PxMala)",ylim = c(-.6, .1))

### Barker
acf_isb <- acf(output_chain[[3]][,1], plot = FALSE)$acf
acf_pxb <- acf(output_chain[[4]][,1], plot = FALSE)$acf

plot(1:length(acf_isb), acf_isb, col = "blue", type = "l",
     xlab = "Lag", ylab = "Autocorrelation", ylim = c(-0.2, 1))
lines(1:length(acf_pxb), acf_pxb, col = "red", type = "l")
legend("bottomright", c("MYBarker", "PxBarker"), lty = 1,
       col = c("blue", "red"), cex = 0.75, bty = "n")

for (i in 2:dim) 
{
  acf_isb <- acf(output_chain[[3]][,i], plot = FALSE)$acf
  acf_pxb <- acf(output_chain[[4]][,i], plot = FALSE)$acf
  lines(1:length(acf_isb), acf_isb, col = "blue", type = "l")
  lines(1:length(acf_pxb), acf_pxb, col = "red", type = "l")
}

### HMC
acf_is_hmc <- acf(output_chain[[5]][,1], plot = FALSE)$acf
acf_pxhmc <- acf(output_chain[[6]][,1], plot = FALSE)$acf

plot(1:length(acf_is_hmc), acf_is_hmc, col = "blue", type = "l",
     xlab = "Lag", ylab = "Autocorrelation", ylim = c(-0.2, 1))
lines(1:length(acf_pxhmc), acf_pxhmc, col = "red", type = "l")
legend("bottomright", c("MYHMC", "PxHMC"), lty = 1,
       col = c("blue", "red"), cex = 0.75, bty = "n")
for (i in 2:dim) 
{
  acf_is_hmc <- acf(output_chain[[3]][,i], plot = FALSE)$acf
  acf_pxhmc <- acf(output_chain[[4]][,i], plot = FALSE)$acf
  lines(1:length(acf_is_hmc), acf_is_hmc, col = "blue", type = "l")
  lines(1:length(acf_pxhmc), acf_pxhmc, col = "red", type = "l")
}
dev.off()

##################  Boxplots of marginal efficiency ##################

#### Marginal variance comparison for IS vs PxMALA
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
rel_eff_mat_mala <- margvar_pxm/margvar_ism

#### Marginal variance comparison for IS vs PxBarker
#x <- c(1:100)
#dim <- 51 
margvar_isb <- matrix(0, nrow = 100, ncol = dim)
margvar_pxb <- matrix(0, nrow = 100, ncol = dim)
for (i in 1:100)
{
  asympmat_isb <- matrix(unlist(output_poisson[[i]][[3]]), nrow = dim, ncol = dim, byrow = T)
  asympmat_pxb <- matrix(unlist(output_poisson[[i]][[4]]), nrow = dim, ncol = dim, byrow = T)
  for (j in 1:dim)
  {
    margvar_isb[i,j] <- asympmat_isb[j,j]
    margvar_pxb[i,j] <- asympmat_pxb[j,j]
  }
}
rel_eff_mat_bark <- margvar_pxb/margvar_isb

#### Marginal variance comparison for IS vs PxHMC

margvar_ish <- matrix(0, nrow = 100, ncol = dim)
margvar_pxh <- matrix(0, nrow = 100, ncol = dim)
for (i in 1:100)
{
  asympmat_ish <- matrix(unlist(output_poisson[[i]][5]), nrow = dim, ncol = dim, byrow = T)
  asympmat_pxh <- matrix(unlist(output_poisson[[i]][6]), nrow = dim, ncol = dim, byrow = T)
  for (j in 1:dim)
  {
    margvar_ish[i,j] <- asympmat_ish[j,j]
    margvar_pxh[i,j] <- asympmat_pxh[j,j]
  }
}
rel_eff_mat_hmc <- margvar_pxh/margvar_ish

pdf(file = "plots/boxplot_mala.pdf", width = 12, height = 8)
boxplot(rel_eff_mat_mala, use.cols = TRUE, xlab = "Coordinate", 
        ylab = "Average relative efficiency")
dev.off()
pdf(file = "plots/boxplot_bark.pdf", width = 12, height = 8)
boxplot(rel_eff_mat_bark, use.cols = TRUE, xlab = "Coordinate", 
        ylab = "Average relative efficiency")
dev.off()
pdf(file = "plots/boxplot_hmc.pdf", width = 12, height = 8)
boxplot(rel_eff_mat_hmc, use.cols = TRUE, xlab = "Coordinate", 
        ylab = "Average relative efficiency")
dev.off()

###############  Histogram  ###############

avg_rel_eff_mala <- apply(rel_eff_mat_mala, 2, mean)
avg_rel_eff_bark <- apply(rel_eff_mat_bark, 2, mean)
avg_rel_eff_hmc <- apply(rel_eff_mat_hmc, 2, mean)

pdf(file = "plots/hist_pois_mala.pdf", width = 10, height = 8)
hist(avg_rel_eff_mala, breaks = 10, main = NULL, xlab = "Average relative efficiency")
dev.off()
pdf(file = "plots/hist_pois_bark.pdf", width = 10, height = 8)
hist(avg_rel_eff_bark, breaks = 10, main = NULL, xlab = "Average relative efficiency")
dev.off()
pdf(file = "plots/hist_pois_hmc.pdf", width = 10, height = 8)
hist(avg_rel_eff_hmc, breaks = 10, main = NULL, xlab = "Average relative efficiency")
dev.off()
