################################################################################
################## Trendfiltering example output visualisation #################
################################################################################
rm(list = ls())
library(ggplot2)
source("TF_functions.R")
load("single_chain_mala.Rdata")
load("single_chain_hmc.Rdata")
load("single_chain_log_weights.Rdata")

## weights
mala_wts <- exp(log_wts[[1]])
mean(mala_wts)^2/mean(mala_wts^2)

hmc_wts <- exp(log_wts[[2]])
mean(hmc_wts)^2/mean(hmc_wts^2)

################ ACF plots ################
dim <- 100


### MALA #######

pdf("plots/tf_acf_MALA.pdf", height = 6, width = 6)
par(mfrow = c(1,1))
lag.max <- 30
acf_ism <- acf(output_single_run_mala[[1]][,1], plot = FALSE, lag.max = lag.max)$acf
acf_pxm <- acf(output_single_run_mala[[2]][,1], plot = FALSE, lag.max = lag.max)$acf

diff.acf <- matrix(0, ncol = dim, nrow = lag.max + 1)
diff.acf[,1] <- acf_ism - acf_pxm

for (i in 2:100) 
{
  acf_ism <- acf(output_single_run_mala[[1]][,i], plot = FALSE, lag.max = lag.max)$acf
  acf_pxm <- acf(output_single_run_mala[[2]][,i], plot = FALSE, lag.max = lag.max)$acf
  diff.acf[,i] <- acf_ism - acf_pxm
  # lines(1:length(acf_ism), acf_ism, col = "blue", type = "l")
  # lines(1:length(acf_pxm), acf_pxm, col = "red", type = "l")
}

# Make boxplot of ACFs
boxplot(t(diff.acf),
        xlab = "Lags", col = "pink",
        ylab = "Difference in ACFs of MALAs",ylim = range(diff.acf),
        names = 0:lag.max, show.names = TRUE, cex.lab = 1.5)
dev.off()

##### HMC  ##### 

pdf("plots/tf_acf_HMC.pdf", height = 6, width = 6)
par(mfrow = c(1,1))
acf_is_hmc <- acf(output_single_run_hmc[[1]][,1], plot = FALSE, lag.max = lag.max)$acf
acf_pxhmc <- acf(output_single_run_hmc[[2]][,1], plot = FALSE, lag.max = lag.max)$acf

diff.acf <- matrix(0, ncol = dim, nrow = lag.max + 1)
diff.acf[,1] <- acf_is_hmc - acf_pxhmc

for (i in 2:100) 
{
  acf_is_hmc <- acf(output_single_run_hmc[[1]][,i], plot = FALSE, lag.max = lag.max)$acf
  acf_pxhmc <- acf(output_single_run_hmc[[2]][,i], plot = FALSE, lag.max = lag.max)$acf
  diff.acf[,i] <- acf_is_hmc - acf_pxhmc
}

# Make boxplot of ACFs
boxplot(t(diff.acf),
        xlab = "Lags", col = "pink",
        ylab = "Difference in ACF of HMCs",ylim = range(diff.acf),
        names = 0:lag.max, show.names = TRUE, cex.lab = 1.5)

dev.off()



################## Quantiles for pi ##################

level <- 0.025
### MALA
upper_quant_mala_pi <- quantile_func(output_single_run_mala[[1]], exp(log_wts[[1]]), level)[[1]]
lower_quant_mala_pi <- quantile_func(output_single_run_mala[[1]], exp(log_wts[[1]]), level)[[2]]
post_mean_mala_pi <- post_mean_fn(output_single_run_mala[[1]], log_wts[[1]])

### HMC
upper_quant_hmc_pi <- quantile_func(output_single_run_hmc[[1]], exp(log_wts[[2]]), level)[[1]]
lower_quant_hmc_pi <- quantile_func(output_single_run_hmc[[1]], exp(log_wts[[2]]), level)[[2]]
post_mean_hmc_pi <- post_mean_fn(output_single_run_hmc[[1]], log_wts[[2]])

######################## Quantiles for pi^lambda ########################

### MALA
upper_quant_mala_pilambda <- numeric(length = length(y))
lower_quant_mala_pilambda <- numeric(length = length(y))
post_mean_mala_pilambda <- colMeans(output_single_run_mala[[1]])
for(i in 1:length(y))
{
  upper_quant_mala_pilambda[i] <- quantile(output_single_run_mala[[1]][,i], probs = 0.975)
  lower_quant_mala_pilambda[i] <- quantile(output_single_run_mala[[1]][,i], probs = 0.025)
}

### HMC
upper_quant_hmc_pilambda <- numeric(length = length(y))
lower_quant_hmc_pilambda <- numeric(length = length(y))
post_mean_hmc_pilambda <- colMeans(output_single_run_hmc[[1]])
for(i in 1:length(y))
{
  upper_quant_hmc_pilambda[i] <- quantile(output_single_run_hmc[[1]][,i], probs = 0.975)
  lower_quant_hmc_pilambda[i] <- quantile(output_single_run_hmc[[1]][,i], probs = 0.025)
}

pdf(file = "plots/tf_quantiles_is_mala.pdf", width = 12, height = 6)
dataset <- data.frame(x, y, lower_quant_mala_pi, upper_quant_mala_pi, post_mean_mala_pi)

plot <- ggplot(dataset, aes(x = x)) + 
  geom_ribbon(aes(ymin = lower_quant_mala_pi, ymax = upper_quant_mala_pi,
                  fill = "95% Credible Interval"), alpha = 0.3) +
  geom_point(aes(y = y, color = "Observed")) +
  geom_line(aes(y = post_mean_mala_pi, color = "Posterior Mean"), size = 1.2) +
  scale_color_discrete(name = NULL) +
  scale_fill_discrete(name = NULL) +
  labs(x = "Index", y = "y") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "grey90", color = NA),
    plot.background = element_rect(fill = "grey90", color = NA),
    
    # Axis labels: "Index" and "y"
    axis.title = element_text(size = 16),
    
    # Axis tick labels
    axis.text = element_text(size = 13),
    
    # Legend text
    legend.text = element_text(size = 14),
    
    # Legend position
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black")
  )

print(plot)
dev.off()







##################  Replication study ##################
load("output_mala.Rdata")
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


###### Component-wise average relative efficiency  ######
mu_mala <- colMeans(rel_eff_mat_mala, na.rm = TRUE)
lower_mala <- mu_mala - 1.96*apply(rel_eff_mat_mala, 2, sd)/sqrt(100)
upper_mala <- mu_mala + 1.96*apply(rel_eff_mat_mala, 2, sd)/sqrt(100)


pdf(file = "plots/tf_unc_mala.pdf", width = 10, height = 5)
index <- 1:ncol(rel_eff_mat_mala)
plot(index, mu_mala, pch = 16, col = "pink", ylim = c(0.5, 4),
     xlab = "Component", ylab = "Average efficiency", cex.lab = 1.5)
segments(index, lower_mala, index, mu_mala)
segments(index, upper_mala, index, mu_mala)
dev.off()

# histogram of average relative efficiencies
pdf("plots/tf_boxeff_mala.pdf", height = 5, width = 6)
boxplot(mu_mala, names = "Relative efficiency for MALA", show.names = TRUE,
        boxwex = .5, col = "pink", horizontal  = TRUE, cex.lab = 1.5)
dev.off()

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


# HMC
mu_hmc <- colMeans(rel_eff_mat_hmc, na.rm = TRUE)
lower_hmc <- mu_hmc - 1.96*apply(rel_eff_mat_hmc, 2, sd)/sqrt(100)
upper_hmc <- mu_hmc + 1.96*apply(rel_eff_mat_hmc, 2, sd)/sqrt(100)


###### Component-wise average relative efficiency  ######
pdf(file = "plots/tf_unc_hmc.pdf", width = 10, height = 5)
index <- 1:ncol(rel_eff_mat_mala)
plot(index, mu_hmc, pch = 16, col = "pink", ylim = c(0.5, 40),
     xlab = "Component", ylab = "Average efficiency", cex.lab = 1.5)
segments(index, lower_hmc, index, mu_hmc)
segments(index, upper_hmc, index, mu_hmc)

dev.off() 

## Histrogram of average relative efficiencies
pdf("plots/tf_boxeff_hmc.pdf", height = 5, width = 6)
boxplot(mu_hmc, names = "Relative efficiency for HMC", show.names = TRUE,
        boxwex = .5, col = "pink", horizontal = TRUE, cex.lab = 1.5)
dev.off()




