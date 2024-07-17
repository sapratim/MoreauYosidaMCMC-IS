
source("nuclear_norm_functions.R")
load("output_single chain.Rdata")
load("output_nucl_norm.Rdata")

subset <- 100
rand <- sample(c(1:length(y)), subset)
################ ACF plots ################

#### MALA

pdf("plots/acf_nnorm.pdf", height = 8, width = 15)
par(mfrow = c(1,3))

acf_ism <- acf(output_single[[1]][,rand[1]], plot = FALSE)$acf
acf_pxm <- acf(output_single[[2]][,rand[1]], plot = FALSE)$acf

plot(1:length(acf_ism), acf_ism, col = "blue", type = "l",
     xlab = "Lag", ylab = "Autocorrelation", ylim = c(-0.2, 1), main = "MALA")
lines(1:length(acf_pxm), acf_pxm, col = "red", type = "l")
legend("bottomright", c("MYMALA", "PxMALA"), lty = 1,
       col = c("blue", "red"), cex = 0.75, bty = "n")

for (i in 2:100) 
{
  acf_ism <- acf(output_single[[1]][,rand[i]], plot = FALSE)$acf
  acf_pxm <- acf(output_single[[2]][,rand[i]], plot = FALSE)$acf
  lines(1:length(acf_ism), acf_ism, col = "blue", type = "l")
  lines(1:length(acf_pxm), acf_pxm, col = "red", type = "l")
}

#### Barker

acf_isb <- acf(output_single[[3]][,rand[1]], plot = FALSE)$acf
acf_pxb <- acf(output_single[[4]][,rand[1]], plot = FALSE)$acf

plot(1:length(acf_isb), acf_isb, col = "blue", type = "l",
     xlab = "Lag", ylab = "Autocorrelation", ylim = c(-0.2, 1), main = "Barker")
lines(1:length(acf_pxb), acf_pxb, col = "red", type = "l")
legend("bottomright", c("MYMALA", "PxMALA"), lty = 1,
       col = c("blue", "red"), cex = 0.75, bty = "n")

for (i in 2:100) 
{
  acf_isb <- acf(output_single[[3]][,rand[i]], plot = FALSE)$acf
  acf_pxb <- acf(output_single[[4]][,rand[i]], plot = FALSE)$acf
  lines(1:length(acf_isb), acf_isb, col = "blue", type = "l")
  lines(1:length(acf_pxb), acf_pxb, col = "red", type = "l")
}

####  HMC

acf_is_hmc <- acf(output_single[[5]][,rand[1]], plot = FALSE)$acf
acf_pxhmc <- acf(output[[6]][,rand[1]], plot = FALSE)$acf

plot(1:length(acf_is_hmc), acf_is_hmc, col = "blue", type = "l",
     xlab = "Lag", ylab = "Autocorrelation", ylim = c(-0.2, 1), main = "HMC")
lines(1:length(acf_pxhmc), acf_pxhmc, col = "red", type = "l")
legend("bottomright", c("MYHMC", "PxHMC"), lty = 1,
       col = c("blue", "red"), cex = 0.75, bty = "n")
for (i in 2:100) 
{
  acf_is_hmc <- acf(output_single[[5]][,rand[i]], plot = FALSE)$acf
  acf_pxhmc <- acf(output_single[[6]][,rand[i]], plot = FALSE)$acf
  lines(1:length(acf_is_hmc), acf_is_hmc, col = "blue", type = "l")
  lines(1:length(acf_pxhmc), acf_pxhmc, col = "red", type = "l")
}
dev.off()

################  Boxplots visualisation  ################

mar_eff_mala <- matrix(0, nrow = subset, ncol = length(y))
mar_eff_bark <- matrix(0, nrow = subset, ncol = length(y))
mar_eff_hmc <- matrix(0, nrow = subset, ncol = length(y))

for (i in 1:subset) {
  mar_eff_mala[i,] <- as.numeric(unlist(output[[i]][2]))/as.numeric(unlist(output[[i]][1]))
  mar_eff_bark[i,] <- as.numeric(unlist(output[[i]][4]))/as.numeric(unlist(output[[i]][3]))
  mar_eff_hmc[i,] <- as.numeric(unlist(output[[i]][6]))/as.numeric(unlist(output[[i]][5]))
}

pdf(file = "plots/boxplot_nnorm_mala.pdf", width = 10, height = 8)
boxplot(mar_eff_mala[,rand], use.cols = TRUE, xlab = "Coordinate", 
        ylab = "Average relative efficiency")
dev.off()
pdf(file = "plots/boxplot_nnorm_bark.pdf", width = 10, height = 8)
boxplot(mar_eff_bark[,rand], use.cols = TRUE, xlab = "Coordinate", 
        ylab = "Average relative efficiency")
dev.off()
pdf(file = "plots/boxplot_nnorm_hmc.pdf", width = 10, height = 8)
boxplot(mar_eff_hmc[,rand], use.cols = TRUE, xlab = "Coordinate", 
        ylab = "Average relative efficiency")
dev.off()

#############  Histogram

avg_rel_eff_mala <- apply(mar_eff_mala, 2, mean)
avg_rel_eff_bark <- apply(mar_eff_bark, 2, mean)
avg_rel_eff_hmc <- apply(mar_eff_hmc, 2, mean)

pdf(file = "plots/hist_nnorm_mala.pdf", width = 8, height = 6)
hist(avg_rel_eff_mala, breaks = 30, xlab = "Average relative efficiency", main = NULL)
dev.off()

pdf(file = "plots/hist_nnorm_bark.pdf", width = 8, height = 6)
hist(avg_rel_eff_bark, breaks = 30, xlab = "Average relative efficiency", main = NULL)
dev.off()

pdf(file = "plots/hist_nnorm_hmc.pdf", width = 8, height = 6)
hist(avg_rel_eff_hmc, breaks = 30, xlab = "Average relative efficiency", main = NULL)
dev.off()

