################################################################################
################## Trendfiltering example output visualisation #################
################################################################################
rm(list = ls())
source("TF_functions.R")
load("output_mala.Rdata")
load("output_bark.Rdata")
load("output_hmc.Rdata")
load("single_chain.Rdata")
load("single_chain_log_weights.Rdata")

################ ACF plots ################

pdf("plots/acf_tf.pdf", height = 6, width = 12)
par(mfrow = c(1,3))

### MALA
acf_ism <- acf(output_single_run[[1]][,1], plot = FALSE)$acf
acf_pxm <- acf(output_single_run[[2]][,1], plot = FALSE)$acf

plot(1:length(acf_ism), acf_ism, col = "blue", type = "l",
     xlab = "Lag", ylab = "Autocorrelation", ylim = c(-0.2, 1))
lines(1:length(acf_pxm), acf_pxm, col = "red", type = "l")
legend("bottomright", c("MYMALA", "PxMALA"), lty = 1,
       col = c("blue", "red"), cex = 0.75, bty = "n")

for (i in 2:100) 
{
  acf_ism <- acf(output_single_run[[1]][,i], plot = FALSE)$acf
  acf_pxm <- acf(output_single_run[[2]][,i], plot = FALSE)$acf
  lines(1:length(acf_ism), acf_ism, col = "blue", type = "l")
  lines(1:length(acf_pxm), acf_pxm, col = "red", type = "l")
}

### Barker
acf_isb <- acf(output_single_run[[3]][,1], plot = FALSE)$acf
acf_pxb <- acf(output_single_run[[4]][,1], plot = FALSE)$acf

plot(1:length(acf_isb), acf_isb, col = "blue", type = "l",
     xlab = "Lag", ylab = "Autocorrelation", ylim = c(-0.2, 1))
lines(1:length(acf_pxb), acf_pxb, col = "red", type = "l")
legend("bottomright", c("MYBarker", "PxBarker"), lty = 1,
       col = c("blue", "red"), cex = 0.75, bty = "n")

for (i in 2:100) 
{
  acf_isb <- acf(output_single_run[[3]][,i], plot = FALSE)$acf
  acf_pxb <- acf(output_single_run[[4]][,i], plot = FALSE)$acf
  lines(1:length(acf_isb), acf_isb, col = "blue", type = "l")
  lines(1:length(acf_pxb), acf_pxb, col = "red", type = "l")
}

### HMC
acf_is_hmc <- acf(output_single_run[[5]][,1], plot = FALSE)$acf
acf_pxhmc <- acf(output_single_run[[6]][,1], plot = FALSE)$acf

plot(1:length(acf_is_hmc), acf_is_hmc, col = "blue", type = "l",
     xlab = "Lag", ylab = "Autocorrelation", ylim = c(-0.2, 1))
lines(1:length(acf_pxhmc), acf_pxhmc, col = "red", type = "l")
legend("bottomright", c("MYHMC", "PxHMC"), lty = 1,
       col = c("blue", "red"), cex = 0.75, bty = "n")
for (i in 2:100) 
{
  acf_is_hmc <- acf(output_single_run[[3]][,i], plot = FALSE)$acf
  acf_pxhmc <- acf(output_single_run[[4]][,i], plot = FALSE)$acf
  lines(1:length(acf_is_hmc), acf_is_hmc, col = "blue", type = "l")
  lines(1:length(acf_pxhmc), acf_pxhmc, col = "red", type = "l")
}
dev.off()

##################  Boxplots of marginal efficiency ##################

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
x <- c(1:100)
dim <- 100 
margvar_isb <- matrix(0, nrow = 100, ncol = dim)
margvar_pxb <- matrix(0, nrow = 100, ncol = dim)
for (i in 1:100)
{
  asympmat_isb <- matrix(unlist(output_bark[[i]][[3]]), nrow = dim, ncol = dim, byrow = T)
  asympmat_pxb <- matrix(unlist(output_bark[[i]][[4]]), nrow = dim, ncol = dim, byrow = T)
  for (j in 1:dim)
  {
    margvar_isb[i,j] <- asympmat_ism[j,j]
    margvar_pxb[i,j] <- asympmat_pxm[j,j]
  }
}
rel_eff_mat_bark <- margvar_pxm/margvar_ism

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

###############  Histogram

avg_rel_eff_mala <- apply(rel_eff_mat_mala, 2, mean)
avg_rel_eff_bark <- apply(rel_eff_mat_bark, 2, mean)
avg_rel_eff_hmc <- apply(rel_eff_mat_hmc, 2, mean)

pdf(file = "plots/hist_tf_mala.pdf", width = 10, height = 8)
hist(avg_rel_eff_mala, breaks = 10, main = NULL, xlab = "Average relative efficiency")
dev.off()
pdf(file = "plots/hist_tf_bark.pdf", width = 10, height = 8)
hist(avg_rel_eff_bark, breaks = 10, main = NULL, xlab = "Average relative efficiency")
dev.off()
pdf(file = "plots/hist_tf_hmc.pdf", width = 10, height = 8)
hist(avg_rel_eff_hmc, breaks = 10, main = NULL, xlab = "Average relative efficiency")
dev.off()

################## Quantiles for pi ##################

level <- 0.025
### MALA
upper_quant_mala_pi <- quantile_func(output_single_run[[1]], exp(log_wts[[1]]), level)[[1]]
lower_quant_mala_pi <- quantile_func(output_single_run[[1]], exp(log_wts[[1]]), level)[[2]]
post_med_mala_pi <- quantile_func(output_single_run[[1]], exp(log_wts[[1]]), level)[[3]]

### Barker
upper_quant_bark_pi <- quantile_func(output_single_run[[3]], exp(log_wts[[2]]), level)[[1]]
lower_quant_bark_pi <- quantile_func(output_single_run[[3]], exp(log_wts[[2]]), level)[[2]]
post_med_bark_pi <- quantile_func(output_single_run[[3]], exp(log_wts[[2]]), level)[[3]]

### HMC
upper_quant_hmc_pi <- quantile_func(output_single_run[[5]], exp(log_wts[[3]]), level)[[1]]
lower_quant_hmc_pi <- quantile_func(output_single_run[[5]], exp(log_wts[[3]]), level)[[2]]
post_med_hmc_pi <- quantile_func(output_single_run[[5]], exp(log_wts[[3]]), level)[[3]]

######################## Quantiles for pi^lambda ########################

### MALA
upper_quant_mala_pilambda <- numeric(length = length(y))
lower_quant_mala_pilambda <- numeric(length = length(y))
post_med_mala_pilambda <- numeric(length = length(y))
for(i in 1:length(y))
{
  upper_quant_mala_pilambda[i] <- quantile(output_single_run[[1]][,i], probs = 0.975)
  lower_quant_mala_pilambda[i] <- quantile(output_single_run[[1]][,i], probs = 0.025)
  post_med_mala_pilambda[i] <- quantile(output_single_run[[1]][,i], probs = 0.5)
}

### Barker
upper_quant_bark_pilambda <- numeric(length = length(y))
lower_quant_bark_pilambda <- numeric(length = length(y))
post_med_bark_pilambda <- numeric(length = length(y))
for(i in 1:length(y))
{
  upper_quant_bark_pilambda[i] <- quantile(output_single_run[[3]][,i], probs = 0.975)
  lower_quant_bark_pilambda[i] <- quantile(output_single_run[[3]][,i], probs = 0.025)
  post_med_bark_pilambda[i] <- quantile(output_single_run[[3]][,i], probs = 0.5)
}

### HMC
upper_quant_hmc_pilambda <- numeric(length = length(y))
lower_quant_hmc_pilambda <- numeric(length = length(y))
post_med_hmc_pilambda <- numeric(length = length(y))
for(i in 1:length(y))
{
  upper_quant_hmc_pilambda[i] <- quantile(output_single_run[[5]][,i], probs = 0.975)
  lower_quant_hmc_pilambda[i] <- quantile(output_single_run[[5]][,i], probs = 0.025)
  post_med_hmc_pilambda[i] <- quantile(output_single_run[[5]][,i], probs = 0.5)
}

pdf(file = "plots/tf_quantiles_is_mala.pdf", width = 10, height = 6)
dataset <- data.frame(x, y, lower_quant_mala_pi, upper_quant_mala_pi, post_med_mala_pi)
plot <- ggplot(dataset, aes(x, y,group = )) + geom_point() +
  geom_line(aes(x=c(1:100), y=post_med), col = "red")
conf_bands <- plot + geom_ribbon(aes(ymin = lower_quant_mala_pi, ymax = upper_quant_mala_pi),
                 alpha = 0.3) +labs(x = "variable") + labs(y = "data")
conf_bands
dev.off()
