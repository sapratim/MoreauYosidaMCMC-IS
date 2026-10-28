
source("nuclear_norm_functions.R")

# subset <- 100
rand <- 1:length(y) #sample(c(1:length(y)), subset)
dim <- length(y)

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

# ###################. Boxplots ####################
# 
# boxplot(mar_eff_mala, main = "Column-wise Boxplots", xlab = "Columns",
#         ylab = "Values", col = "lightblue",names = paste0("Column ", 1:ncol(mar_eff_mala)))
# boxplot(mar_eff_hmc, main = "Column-wise Boxplots", xlab = "Columns",
#         ylab = "Values", col = "lightblue",names = paste0("Column ", 1:ncol(mar_eff_hmc)))


# MALA
lower_mala <- apply(mar_eff_mala, 2, quantile, probs = 0.025, na.rm = TRUE)
upper_mala <- apply(mar_eff_mala, 2, quantile, probs = 0.925, na.rm = TRUE)
mu <- colMeans(mar_eff_mala, na.rm = TRUE)

pdf(file = "Rcpp/plots/boxplots_unc_mala.pdf", width = 10, height = 5)
# Plot the band
matplot(1:ncol(mar_eff_mala),cbind(lower_mala, upper_mala, mu),type = "l",lty = c(2, 2, 1), 
        lwd = 2,col = c("blue", "blue", "red"),xlab = "Column",ylab = "Value",
        main = "90% Quantile Band and Mean for MALA")
dev.off()

# HMC
lower_hmc <- apply(mar_eff_hmc, 2, quantile, probs = 0.025, na.rm = TRUE)
upper_hmc <- apply(mar_eff_hmc, 2, quantile, probs = 0.925, na.rm = TRUE)
mu <- colMeans(mar_eff_hmc, na.rm = TRUE)

pdf(file = "Rcpp/plots/boxplots_unc_hmc.pdf", width = 10, height = 5)
# Plot the band
matplot(1:ncol(mar_eff_hmc),cbind(lower_hmc, upper_hmc, mu),type = "l",lty = c(2, 2, 1), 
        lwd = 2,col = c("blue", "blue", "red"),xlab = "Column",ylab = "Value",
        main = "90% Quantile Band and Mean for HMC")
dev.off() 
