################################################################################
################## Trendfiltering example output visualisation #################
################################################################################
rm(list = ls())
library(ggplot2)
source("TF_functions.R")
load("single_chain_mala.Rdata")
#load("single_chain_bark.Rdata")
load("single_chain_hmc.Rdata")
load("single_chain_log_weights.Rdata")

## weights
mala_wts <- exp(log_wts[[1]])
mean(mala_wts)^2/mean(mala_wts^2)

hmc_wts <- exp(log_wts[[3]])
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
        names = 0:lag.max, show.names = TRUE)
dev.off()
######## Barker #################

# acf_isb <- acf(output_single_run_bark[[1]][,1], plot = FALSE, lag.max = lag.max)$acf
# acf_pxb <- acf(output_single_run_bark[[2]][,1], plot = FALSE, lag.max = lag.max)$acf

# diff.acf <- matrix(0, ncol = dim, nrow = lag.max + 1)
# diff.acf[,1] <- acf_isb - acf_pxb

# # plot(1:length(acf_isb), acf_isb, col = "blue", type = "l",
# #      xlab = "Lag", ylab = "Autocorrelation", ylim = c(-0.2, 1))
# # lines(1:length(acf_pxb), acf_pxb, col = "red", type = "l")
# # legend("bottomright", c("MYBarker", "PxBarker"), lty = 1,
# #        col = c("blue", "red"), cex = 0.75, bty = "n")

# for (i in 2:100) 
# {
#   acf_isb <- acf(output_single_run_bark[[1]][,i], plot = FALSE, lag.max = lag.max)$acf
#   acf_pxb <- acf(output_single_run_bark[[2]][,i], plot = FALSE, lag.max = lag.max)$acf
#   diff.acf[,i] <- acf_isb - acf_pxb
#   # lines(1:length(acf_isb), acf_isb, col = "blue", type = "l")
#   # lines(1:length(acf_pxb), acf_pxb, col = "red", type = "l")
# }

# # Make boxplot of ACFs
# boxplot(t(diff.acf),
#         xlab = "Lags", ylab = "Difference in ACF (MYBarker - PxBarker)",ylim = c(-.6, .1))

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
        names = 0:lag.max, show.names = TRUE)

dev.off()



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
rel_eff_mat_mala <- margvar_pxm/margvar_ism

#### Marginal variance comparison for IS vs PxBarker
# x <- c(1:100)
# dim <- 100 
# margvar_isb <- matrix(0, nrow = 100, ncol = dim)
# margvar_pxb <- matrix(0, nrow = 100, ncol = dim)
# for (i in 1:100)
# {
#   asympmat_isb <- matrix(unlist(output_bark[[i]][[3]]), nrow = dim, ncol = dim, byrow = T)
#   asympmat_pxb <- matrix(unlist(output_bark[[i]][[4]]), nrow = dim, ncol = dim, byrow = T)
#   for (j in 1:dim)
#   {
#     margvar_isb[i,j] <- asympmat_isb[j,j]
#     margvar_pxb[i,j] <- asympmat_pxb[j,j]
#   }
# }
# rel_eff_mat_bark <- margvar_pxb/margvar_isb

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
rel_eff_mat_hmc <- margvar_pxh/margvar_ish

# pdf(file = "plots/boxplot_mala.pdf", width = 12, height = 8)
# boxplot(rel_eff_mat_mala, use.cols = TRUE, xlab = "Coordinate", 
#         ylab = "Average relative efficiency")
# dev.off()
# pdf(file = "plots/boxplot_bark.pdf", width = 12, height = 8)
# boxplot(rel_eff_mat_bark, use.cols = TRUE, xlab = "Coordinate", 
#         ylab = "Average relative efficiency")
# dev.off()
# pdf(file = "plots/boxplot_hmc.pdf", width = 12, height = 8)
# boxplot(rel_eff_mat_hmc, use.cols = TRUE, xlab = "Coordinate", 
#         ylab = "Average relative efficiency")
# dev.off()

###############  Histogram

avg_rel_eff_mala <- apply(rel_eff_mat_mala, 2, mean)
# avg_rel_eff_bark <- apply(rel_eff_mat_bark, 2, mean)
avg_rel_eff_hmc <- apply(rel_eff_mat_hmc, 2, mean)

pdf("plots/tf_boxeff_mala.pdf", height = 5, width = 6)
boxplot(avg_rel_eff_mala, names = "Relative efficiency for MALA", show.names = TRUE,
        boxwex = .5, col = "pink", horizontal  = TRUE)
dev.off()

pdf("plots/tf_boxeff_hmc.pdf", height = 5, width = 6)
boxplot(avg_rel_eff_hmc, names = "Relative efficiency for HMC", show.names = TRUE,
        boxwex = .5, col = "pink", horizontal = TRUE)
dev.off()

# pdf(file = "plots/hist_tf_mala.pdf", width = 10, height = 8)
# hist(avg_rel_eff_mala, breaks = 10, main = NULL, xlab = "Average relative efficiency")
# dev.off()
# pdf(file = "plots/hist_tf_bark.pdf", width = 10, height = 8)
# hist(avg_rel_eff_bark, breaks = 10, main = NULL, xlab = "Average relative efficiency")
# dev.off()
# pdf(file = "plots/hist_tf_hmc.pdf", width = 10, height = 8)
# hist(avg_rel_eff_hmc, breaks = 10, main = NULL, xlab = "Average relative efficiency")
# dev.off()

################## Quantiles for pi ##################

level <- 0.025
### MALA
upper_quant_mala_pi <- quantile_func(output_single_run_mala[[1]], exp(log_wts[[1]]), level)[[1]]
lower_quant_mala_pi <- quantile_func(output_single_run_mala[[1]], exp(log_wts[[1]]), level)[[2]]
post_mean_mala_pi <- post_mean_fn(output_single_run_mala[[1]], log_wts[[1]])

### Barker
# upper_quant_bark_pi <- quantile_func(output_single_run_bark[[1]], exp(log_wts[[2]]), level)[[1]]
# lower_quant_bark_pi <- quantile_func(output_single_run_bark[[1]], exp(log_wts[[2]]), level)[[2]]
# post_med_bark_pi <- quantile_func(output_single_run_bark[[1]], exp(log_wts[[2]]), level)[[3]]

### HMC
upper_quant_hmc_pi <- quantile_func(output_single_run_hmc[[1]], exp(log_wts[[3]]), level)[[1]]
lower_quant_hmc_pi <- quantile_func(output_single_run_hmc[[1]], exp(log_wts[[3]]), level)[[2]]
post_mean_hmc_pi <- post_mean_fn(output_single_run_hmc[[1]], log_wts[[3]])

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

### Barker
# upper_quant_bark_pilambda <- numeric(length = length(y))
# lower_quant_bark_pilambda <- numeric(length = length(y))
# post_med_bark_pilambda <- numeric(length = length(y))
# for(i in 1:length(y))
# {
#   upper_quant_bark_pilambda[i] <- quantile(output_single_run_bark[[1]][,i], probs = 0.975)
#   lower_quant_bark_pilambda[i] <- quantile(output_single_run_bark[[1]][,i], probs = 0.025)
#   post_med_bark_pilambda[i] <- quantile(output_single_run_bark[[1]][,i], probs = 0.5)
# }

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
plot <- ggplot(dataset, aes(x, y,group = )) + geom_point() +
  geom_line(aes(x=c(1:100), y=post_mean_mala_pi), col = "red")
conf_bands <- plot + geom_ribbon(aes(ymin = lower_quant_mala_pi, ymax = upper_quant_mala_pi),
                                 alpha = 0.3) +labs(x = "index") + labs(y = "y")
conf_bands
dev.off()
