
rm(list = ls())
source("nuclear_norm_functions.R")
load("output_dummy.Rdata")
load("output_nucl_norm_final.Rdata")

reps <- 100

avg_mar_eff_func <- function(output_list)
{
  mar_eff_mala <- matrix(0, nrow = reps, ncol = 4096)
  mar_eff_hmc <- matrix(0, nrow = reps, ncol = 4096)
for (i in 1:reps) {
  mar_eff_mala[i,] <- as.numeric(unlist(output_list[[i]][2]))/as.numeric(unlist(output_list[[i]][1]))
  mar_eff_hmc[i,] <- as.numeric(unlist(output_list[[i]][4]))/as.numeric(unlist(output_list[[i]][3]))
  }
  avg_rel_eff_mala <- apply(mar_eff_mala, 2, mean)
  avg_rel_eff_hmc <- apply(mar_eff_hmc, 2, mean)
  avg_eff_list <- list(avg_rel_eff_mala, avg_rel_eff_hmc)
  return(avg_eff_list)
}

################ ACF plots ################

rand <- sample(c(1:500), reps)
pdf("acf_nnorm.pdf", height = 10, width = 6)
par(mfrow = c(2,1))

acf_is <- acf(output[[2]][,rand[1]], plot = FALSE)$acf
acf_pxm <- acf(output[[3]][,rand[1]], plot = FALSE)$acf

plot(1:length(acf_is), acf_is, col = "blue", type = "l",
     xlab = "Lag", ylab = "Autocorrelation", ylim = c(-0.2, 1), main = "MALA")
lines(1:length(acf_pxm), acf_pxm, col = "red", type = "l")
legend("bottomright", c("MYMALA", "PxMALA"), lty = 1,
       col = c("blue", "red"), cex = 0.75, bty = "n")

for (i in 2:100) 
{
  acf_is <- acf(output[[2]][,rand[i]], plot = FALSE)$acf
  acf_pxm <- acf(output[[3]][,rand[i]], plot = FALSE)$acf
  lines(1:length(acf_is), acf_is, col = "blue", type = "l")
  lines(1:length(acf_pxm), acf_pxm, col = "red", type = "l")
}

acf_is_hmc <- acf(output[[4]][,rand[1]], plot = FALSE)$acf
acf_pxhmc <- acf(output[[5]][,rand[1]], plot = FALSE)$acf

plot(1:length(acf_is_hmc), acf_is_hmc, col = "blue", type = "l",
     xlab = "Lag", ylab = "Autocorrelation", ylim = c(-0.2, 1), main = "HMC")
lines(1:length(acf_pxhmc), acf_pxhmc, col = "red", type = "l")
legend("bottomright", c("MYHMC", "PxHMC"), lty = 1,
       col = c("blue", "red"), cex = 0.75, bty = "n")
for (i in 2:100) 
{
  acf_is_hmc <- acf(output[[4]][,rand[i]], plot = FALSE)$acf
  acf_pxhmc <- acf(output[[5]][,rand[i]], plot = FALSE)$acf
  lines(1:length(acf_is_hmc), acf_is_hmc, col = "blue", type = "l")
  lines(1:length(acf_pxhmc), acf_pxhmc, col = "red", type = "l")
}
dev.off()

#############  Histogram

freq_data <- avg_mar_eff_func(output_final)
pdf(file = "hist_nnorm_mala.pdf", width = 8, height = 6)
hist(freq_data[[1]], breaks = 30, xlab = "Average relative efficiency", main = NULL)
dev.off()
pdf(file = "hist_nnorm_hmc.pdf", width = 8, height = 6)
hist(freq_data[[2]], breaks = 30, xlab = "Average relative efficiency", main = NULL)
dev.off()

################  Boxplots visualisation  ################

mar_eff_mala <- matrix(0, nrow = reps, ncol = 4096)
mar_eff_hmc <- matrix(0, nrow = reps, ncol = 4096)

for (i in 1:reps) {
  mar_eff_mala[i,] <- as.numeric(unlist(output_final[[i]][2]))/as.numeric(unlist(output_final[[i]][1]))
  mar_eff_hmc[i,] <- as.numeric(unlist(output_final[[i]][4]))/as.numeric(unlist(output_final[[i]][3]))
  }

pdf(file = "boxplot_nnorm_mala.pdf", width = 10, height = 8)
boxplot(mar_eff_mala[,rand], use.cols = TRUE, xlab = "Coordinate", 
        ylab = "Average relative efficiency")
dev.off()
pdf(file = "boxplot_nnorm_hmc.pdf", width = 10, height = 8)
boxplot(mar_eff_hmc[,rand], use.cols = TRUE, xlab = "Coordinate", 
        ylab = "Average relative efficiency")
dev.off()
