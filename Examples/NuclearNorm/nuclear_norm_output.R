
################################################################################
################## Nuclear Norm example output visualisation #################
################################################################################

source("nuclear_norm_functions.R")

# Files below are too big to be on GitHub. Please run
# source("single_chain.R")
# to generate the files. Run time roughly 10 hours
load("output_single_chain_mala.Rdata")
load("output_single_chain_hmc.Rdata")


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
load("output_nucl_norm.Rdata")


mar_eff_mala <- matrix(0, nrow = 100, ncol = length(y))
mar_eff_hmc <- matrix(0, nrow = 100, ncol = length(y))

for (i in 1:100) {
  mar_eff_mala[i,] <- as.numeric(unlist(output[[i]][2]))/as.numeric(unlist(output[[i]][1]))
  mar_eff_hmc[i,] <- as.numeric(unlist(output[[i]][6]))/as.numeric(unlist(output[[i]][5]))
}



#############  Histogram

avg_rel_eff_mala <- apply(mar_eff_mala, 2, mean)
# avg_rel_eff_bark <- apply(mar_eff_bark, 2, mean)
avg_rel_eff_hmc <- apply(mar_eff_hmc, 2, mean)

avg_eff <- cbind(avg_rel_eff_mala, avg_rel_eff_hmc)
colnames(avg_eff) <- c("MALA", "HMC")


pdf("plots/nn_boxeff.pdf", height = 5, width = 8)
boxplot(avg_eff, ylab = "Relative efficiency", xaxt = "n",
  boxwex = .5, col = "pink", horizontal  = TRUE, ylim = c(1,2.5))
axis(1, at = seq(1, 3, by = .5))
dev.off()

