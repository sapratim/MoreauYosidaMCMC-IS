################################################################################
################## Trendfiltering example output visualisation #################
################################################################################
rm(list = ls())
source("IS_trendf_functions_Pereyra.R")
load("output_mala.Rdata")
load("output_hmc.Rdata")
load("single_chain.Rdata")
load("single_chain_log_weights.Rdata")

################ ACF plots ################

pdf("acf_tf.pdf", height = 6, width = 10)
par(mfrow = c(1,2))

acf_is <- acf(output_single_run[[1]][,1], plot = FALSE)$acf
acf_pxm <- acf(output_single_run[[2]][,1], plot = FALSE)$acf

plot(1:length(acf_is), acf_is, col = "blue", type = "l",
     xlab = "Lag", ylab = "Autocorrelation", ylim = c(-0.2, 1))
lines(1:length(acf_pxm), acf_pxm, col = "red", type = "l")
legend("bottomright", c("MYMALA", "PxMALA"), lty = 1,
       col = c("blue", "red"), cex = 0.75, bty = "n")

for (i in 2:100) 
{
  acf_is <- acf(output_single_run[[1]][,i], plot = FALSE)$acf
  acf_pxm <- acf(output_single_run[[2]][,i], plot = FALSE)$acf
  lines(1:length(acf_is), acf_is, col = "blue", type = "l")
  lines(1:length(acf_pxm), acf_pxm, col = "red", type = "l")
}

acf_is_hmc <- acf(output_single_run[[3]][,1], plot = FALSE)$acf
acf_pxhmc <- acf(output_single_run[[4]][,1], plot = FALSE)$acf

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
  asympmat_ism <- matrix(unlist(output[[i]][[4]]), nrow = dim, ncol = dim, byrow = T)
  asympmat_pxm <- matrix(unlist(output[[i]][[5]]), nrow = dim, ncol = dim, byrow = T)
  for (j in 1:dim)
  {
    margvar_ism[i,j] <- asympmat_ism[j,j]
    margvar_pxm[i,j] <- asympmat_pxm[j,j]
  }
}

rel_eff_mat_mala <- margvar_pxm/margvar_ism

#### Marginal variance comparison for IS vs PxHMC
dim <- 100 
margvar_ish <- matrix(0, nrow = 100, ncol = dim)
margvar_pxh <- matrix(0, nrow = 100, ncol = dim)
for (i in 1:100)
{
  asympmat_ish <- matrix(unlist(output_hmc[[i]][4]), nrow = dim, ncol = dim, byrow = T)
  asympmat_pxh <- matrix(unlist(output_hmc[[i]][5]), nrow = dim, ncol = dim, byrow = T)
  for (j in 1:dim)
  {
    margvar_ish[i,j] <- asympmat_ish[j,j]
    margvar_pxh[i,j] <- asympmat_pxh[j,j]
  }
}

rel_eff_mat_hmc <- margvar_pxh/margvar_ish

pdf(file = "boxplot_tf.pdf", width = 12, height = 8)
boxplot(rel_eff_mat_mala, use.cols = TRUE, xlab = "Coordinate", 
        ylab = "Average relative efficiency")
boxplot(rel_eff_mat_hmc, use.cols = TRUE, xlab = "Coordinate", 
        ylab = "Average relative efficiency")
dev.off()

###############  Histogram

avg_rel_eff_mala <- apply(rel_eff_mat_mala, 2, mean)
avg_rel_eff_hmc <- apply(rel_eff_mat_hmc, 2, mean)

pdf(file = "hist_tf_mala.pdf", width = 10, height = 8)
hist(avg_rel_eff_mala, breaks = 10, main = NULL, xlab = "Average relative efficiency")
dev.off()
pdf(file = "hist_tf_hmc.pdf", width = 10, height = 8)
hist(avg_rel_eff_hmc, breaks = 10, main = NULL, xlab = "Average relative efficiency")
dev.off()

################## Quantiles for pi ##################

### MALA
augm_mat <- cbind(output_single_run[[1]],log_wts[[1]])

upper_quant_mala_pi <- numeric(length = length(y))
lower_quant_mala_pi <- numeric(length = length(y))
post_med_mala_pi <- numeric(length = length(y))
signif_level <- 0.025

for (i in 1:length(y)) 
{
  initial_mat <- quant(i, augm_mat)
  true_wts <- exp(initial_mat[,2])
  #   mat_sum <- apply(initial_mat, 2, sum)
  wts_prop <- true_wts/sum(true_wts)
  final_mat <- cbind(initial_mat[1,], cumsum(wts_prop))
  lower_index <- min(which(final_mat[,2] >= signif_level))
  upper_index <- min(which(final_mat[,2] >= (1 - signif_level)))
  med_index <- min(which(final_mat[,2] >= 0.5))
  upper_quant_mala_pi[i] <- initial_mat[upper_index,1]
  lower_quant_mala_pi[i] <- initial_mat[lower_index,1]
  post_med_mala_pi[i] <- initial_mat[med_index,1]
}

### HMC
augm_mat <- cbind(output_single_run[[3]],log_wts[[2]])

upper_quant_hmc_pi <- numeric(length = length(y))
lower_quant_hmc_pi <- numeric(length = length(y))
post_med_hmc_pi <- numeric(length = length(y))
signif_level <- 0.025

for (i in 1:length(y)) 
{
  initial_mat <- quant(i, augm_mat)
  true_wts <- exp(initial_mat[,2])
  #   mat_sum <- apply(initial_mat, 2, sum)
  wts_prop <- true_wts/sum(true_wts)
  final_mat <- cbind(initial_mat[1,], cumsum(wts_prop))
  lower_index <- min(which(final_mat[,2] >= signif_level))
  upper_index <- min(which(final_mat[,2] >= (1 - signif_level)))
  med_index <- min(which(final_mat[,2] >= 0.5))
  upper_quant_hmc_pi[i] <- initial_mat[upper_index,1]
  lower_quant_hmc_pi[i] <- initial_mat[lower_index,1]
  post_med_hmc_pi[i] <- initial_mat[med_index,1]
}

######################## Quantiles for pi^lambda ########################

### MALA
upper_quant_pi_lambda <- numeric(length = length(y))
lower_quant_pi_lambda <- numeric(length = length(y))
post_med_pi_lambda <- numeric(length = length(y))
for(i in 1:length(y))
{
  upper_quant_pi_lambda[i] <- quantile(output_single_run[[1]][,i], probs = 0.975)
  lower_quant_pi_lambda[i] <- quantile(output_single_run[[1]][,i], probs = 0.025)
  post_med_pi_lambda[i] <- quantile(output_single_run[[1]][,i], probs = 0.5)
}

### HMC
upper_quant_pi_lambda.hmc <- numeric(length = length(y))
lower_quant_pi_lambda.hmc <- numeric(length = length(y))
post_med_pi_lambda.hmc <- numeric(length = length(y))
for(i in 1:length(y))
{
  upper_quant_pi_lambda.hmc[i] <- quantile(output_single_run[[3]][,i], probs = 0.975)
  lower_quant_pi_lambda.hmc[i] <- quantile(output_single_run[[3]][,i], probs = 0.025)
  post_med_pi_lambda.hmc[i] <- quantile(output_single_run[[3]][,i], probs = 0.5)
}

#######  Quantiles for pi using Px
## MALA
upper_quant_pi_lambda_pxm <- numeric(length = length(y))
lower_quant_pi_lambda_pxm <- numeric(length = length(y))
post_med_pi_lambda_pxm <- numeric(length = length(y))
for(i in 1:length(y))
{
  upper_quant_pi_lambda_pxm[i] <- quantile(output_single_run[[2]][,i], probs = 0.975)
  lower_quant_pi_lambda_pxm[i] <- quantile(output_single_run[[2]][,i], probs = 0.025)
  post_med_pi_lambda_pxm[i] <- quantile(output_single_run[[2]][,i], probs = 0.5)
}

# pdf(file = "tf_quantiles_mala.pdf", width = 12, height = 6)
# par(mfrow = c(1,2))
# 
# ## pi
# post_med <- post_med_mala_pi
# upper_quant <- upper_quant_mala_pi
# lower_quant <- lower_quant_mala_pi
# dataset <- data.frame(x, y, lower_quant, upper_quant, post_med)
# plot <- ggplot(dataset, aes(x, y,group = )) + geom_point() +
#   geom_line(aes(x=c(1:100), y=post_med), col = "red")
# conf_bands <- plot + geom_ribbon(aes(ymin = lower_quant, ymax = upper_quant), alpha = 0.3) +
#   ggtitle("Piecewise linear model") + labs(x = "variable") + labs(y = "data")
# conf_bands
# 
# ## pi lambda
# post_med <- post_med_pi_lambda
# upper_quant <- upper_quant_pi_lambda
# lower_quant <- lower_quant_pi_lambda
# dataset <- data.frame(x, y, lower_quant, upper_quant, post_med)
# plot <- ggplot(dataset, aes(x, y,group = )) + geom_point() +
#   geom_line(aes(x=c(1:100), y=post_med), col = "red")
# conf_bands <- plot + geom_ribbon(aes(ymin = lower_quant, ymax = upper_quant), alpha = 0.3) +
#   ggtitle("Piecewise linear model (pi^lambda)") + labs(x = "variable") + labs(y = "data")
# conf_bands
# dev.off()
# 
# 
# 
# pdf(file = "tf_quantiles_hmc.pdf", width = 12, height = 6)
# par(mfrow = c(1,2))
# 
# ## pi
# post_med <- post_med_hmc_pi
# upper_quant <- upper_quant_hmc_pi
# lower_quant <- lower_quant_hmc_pi
# dataset <- data.frame(x, y, lower_quant, upper_quant, post_med)
# plot <- ggplot(dataset, aes(x, y,group = )) + geom_point() +
#   geom_line(aes(x=c(1:100), y=post_med), col = "red")
# conf_bands <- plot + geom_ribbon(aes(ymin = lower_quant, ymax = upper_quant), alpha = 0.3) +
#   ggtitle("Piecewise linear model(pi)") + labs(x = "variable") + labs(y = "data")
# conf_bands
# 
# ## pi lambda
# post_med <- post_med_pi_lambda.hmc
# upper_quant <- upper_quant_pi_lambda.hmc
# lower_quant <- lower_quant_pi_lambda.hmc
# dataset <- data.frame(x, y, lower_quant, upper_quant, post_med)
# plot <- ggplot(dataset, aes(x, y,group = )) + geom_point() +
#   geom_line(aes(x=c(1:100), y=post_med), col = "red")
# conf_bands <- plot + geom_ribbon(aes(ymin = lower_quant, ymax = upper_quant), alpha = 0.3) +
#   ggtitle("Piecewise linear model (pi^lambda)") + labs(x = "variable") + labs(y = "data")
# conf_bands
# dev.off()
# 

pdf(file = "tf_quantiles_is_mala.pdf", width = 10, height = 6)

## pi
post_med <- post_med_mala_pi
upper_quant <- upper_quant_mala_pi
lower_quant <- lower_quant_mala_pi
dataset <- data.frame(x, y, lower_quant, upper_quant, post_med)
plot <- ggplot(dataset, aes(x, y,group = )) + geom_point() +
  geom_line(aes(x=c(1:100), y=post_med), col = "red")
conf_bands <- plot + geom_ribbon(aes(ymin = lower_quant, ymax = upper_quant), alpha = 0.3) +
    labs(x = "variable") + labs(y = "data")
conf_bands
dev.off()

pdf(file = "tf_quantiles_px_mala.pdf", width = 10, height = 6)
## Px for pi
post_med <- post_med_pi_lambda_pxm
upper_quant <- upper_quant_pi_lambda_pxm
lower_quant <- lower_quant_pi_lambda_pxm
dataset <- data.frame(x, y, lower_quant, upper_quant, post_med)
plot <- ggplot(dataset, aes(x, y,group = )) + geom_point() +
  geom_line(aes(x=c(1:100), y=post_med), col = "red")
conf_bands <- plot + geom_ribbon(aes(ymin = lower_quant, ymax = upper_quant), alpha = 0.3) +
      labs(x = "variable") + labs(y = "data")
conf_bands
dev.off()
