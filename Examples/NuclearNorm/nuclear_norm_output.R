
################################################################################
################## Nuclear Norm example output visualisation #################
################################################################################

# Run from the project root.
source("Rcpp/nuclear_norm_functions.R")
y <- generate_checkerboard(seed = 8024248)$y

# Files below are too big to be on GitHub. Please run
# source("Rcpp/single_chain.R")
# to generate the files.
load("Rcpp/output_single_chain_mala.Rdata")
load("Rcpp/output_single_chain_hmc.Rdata")


# subset <- 100
rand <- 1:length(y) #sample(c(1:length(y)), subset)
dim <- length(y)
################ ACF plots ################

#### MALA
pdf("plots/nn_acf_MALA.pdf", height = 6, width = 6)
par(mfrow = c(1,1))

lag.max <- 100
acf_ism <- acf(output_single_mala[[1]][,rand[1]], plot = FALSE, lag.max = lag.max)$acf
acf_pxm <- acf(output_single_mala[[2]][,rand[1]], plot = FALSE, lag.max = lag.max)$acf

diff.acf <- matrix(0, ncol = dim, nrow = lag.max + 1)
diff.acf[,1] <- acf_ism - acf_pxm

for (i in 2:dim) 
{
  if(i %% 1000 == 0) print(i)
  acf_ism <- acf(output_single_mala[[1]][,rand[i]], plot = FALSE, lag.max = lag.max)$acf
  acf_pxm <- acf(output_single_mala[[2]][,rand[i]], plot = FALSE, lag.max = lag.max)$acf
  diff.acf[,i] <- acf_ism - acf_pxm
 }

# Make boxplot of ACFs
boxplot(t(diff.acf),
        xlab = "Lags", col = "pink",
        ylab = "Difference in ACFs of MALAs",ylim = c(-.25, .07),
        names = 0:lag.max, show.names = TRUE, range = 3)
dev.off()


########################  HMC  #############################

pdf("plots/nn_acf_HMC.pdf", height = 6, width = 6)
acf_is_hmc <- acf(output_single_hmc[[1]][,rand[1]], plot = FALSE, lag.max = lag.max)$acf
acf_pxhmc <- acf(output_single_hmc[[2]][,rand[1]], plot = FALSE, lag.max = lag.max)$acf

diff.acf <- matrix(0, ncol = dim, nrow = lag.max + 1)
diff.acf[,1] <- acf_is_hmc - acf_pxhmc

for (i in 2:dim) 
{
  if(i %% 1000 == 0) print(i)
  acf_is_hmc <- acf(output_single_hmc[[1]][,rand[i]], plot = FALSE, lag.max = lag.max)$acf
  acf_pxhmc <- acf(output_single_hmc[[2]][,rand[i]], plot = FALSE, lag.max = lag.max)$acf
  diff.acf[,i] <- acf_is_hmc - acf_pxhmc
}


# Make boxplot of ACFs
boxplot(t(diff.acf),
        xlab = "Lags", col = "pink",
        ylab = "Difference in ACFs of HMCs",ylim = range(diff.acf),
        names = 0:lag.max, show.names = TRUE, range = 3)
dev.off()



################  Boxplots visualisation  ################
####### Replications experiment  #####
load("Rcpp/output_nucl_norm.Rdata")

mar_eff_mala <- matrix(0, nrow = 100, ncol = length(y))
mar_eff_hmc <- matrix(0, nrow = 100, ncol = length(y))

time_ism <- as.numeric(sapply(output, function(x) x[[5]]))
time_pxm <- as.numeric(sapply(output, function(x) x[[6]]))
time_ish <- as.numeric(sapply(output, function(x) x[[7]]))
time_pxh <- as.numeric(sapply(output, function(x) x[[8]]))

for (i in 1:100) {
  mar_eff_mala[i,] <- as.numeric(unlist(output[[i]][2])*time_pxm[i])/
                                     as.numeric(unlist(output[[i]][1])*time_ism[i])
  mar_eff_hmc[i,] <- as.numeric(unlist(output[[i]][4])*time_pxh[i])/
                                     as.numeric(unlist(output[[i]][3])*time_ish[i])
}


# MALA
lower_mala <- apply(mar_eff_mala, 2, quantile, probs = 0.025, na.rm = TRUE)
upper_mala <- apply(mar_eff_mala, 2, quantile, probs = 0.925, na.rm = TRUE)
mu <- colMeans(mar_eff_mala, na.rm = TRUE)

pdf(file = "plots/boxplots_unc_mala.pdf", width = 10, height = 5)
# Plot the band
matplot(1:ncol(mar_eff_mala),cbind(lower_mala, upper_mala, mu),type = "l",lty = c(2, 2, 1), 
        lwd = 2,col = c("blue", "blue", "red"),xlab = "Column",ylab = "Value",
        main = "90% Quantile Band and Mean for MALA")
dev.off()

mu_mala <- colMeans(mar_eff_mala, na.rm = TRUE)
lower_mala <- mu_mala - 1.96*apply(mar_eff_mala, 2, sd)/sqrt(100)
upper_mala <- mu_mala + 1.96*apply(mar_eff_mala, 2, sd)/sqrt(100)


pdf(file = "plots/nn_unc_mala.pdf", width = 10, height = 5)
index <- seq(1, ncol(mar_eff_mala), by = 50)
plot(index, mu_mala[index], pch = 16, col = "pink", ylim = c(1, 1.6),
     xlab = "Component", ylab = "Average efficiency", cex.lab = 1.5)
segments(index, lower_mala[index], index, mu_mala[index])
segments(index, upper_mala[index], index, mu_mala[index])
dev.off()


mu_hmc <- colMeans(mar_eff_hmc, na.rm = TRUE)
lower_hmc <- mu_hmc - 1.96*apply(mar_eff_hmc, 2, sd)/sqrt(100)
upper_hmc <- mu_hmc + 1.96*apply(mar_eff_hmc, 2, sd)/sqrt(100)


pdf(file = "plots/nn_unc_hmc.pdf", width = 10, height = 5)
index <- seq(1, ncol(mar_eff_hmc), by = 50)
plot(index, mu_hmc[index], pch = 16, col = "pink", ylim = c(2, 3),
     xlab = "Component", ylab = "Average efficiency", cex.lab = 1.5)
segments(index, lower_hmc[index], index, mu_hmc[index])
segments(index, upper_hmc[index], index, mu_hmc[index])
dev.off()

#############  Histogram

avg_rel_eff_mala <- apply(mar_eff_mala, 2, mean)
avg_rel_eff_hmc <- apply(mar_eff_hmc, 2, mean)

avg_eff <- cbind(avg_rel_eff_mala, avg_rel_eff_hmc)
colnames(avg_eff) <- c("MALA", "HMC")


pdf("plots/nn_boxeff.pdf", height = 5, width = 8)
boxplot(avg_eff, ylab = "Relative efficiency", xaxt = "n",
  boxwex = .5, col = "pink", horizontal  = TRUE, ylim = c(1,2.8))
axis(1, at = seq(1, 3, by = .5))
dev.off()

