source("IS_trendf_functions_Pereyra.R")
load("single_chain.Rdata")

lamb_coeff <- 0.001
D_mat <- getD(k=1, n=1e2, x)   #  D matrix
alpha_hat <- 5   # obtained from the first dataset
sigma2_hat <- 9  # obtained from the first dataset
k <- 1

wts_func <- function(beta,lambda,alpha,sigma2,k,grid)
{
  g_val <- alpha*sum(abs(D_mat%*%beta)) + sum((y - beta)^2)/(2*sigma2)
  proxval <- prox_func(beta, lambda, alpha, sigma2, k, grid)
  g_lambda_val <- prox_arg(proxval, beta, lambda=lambda, y, sigma2, alpha)
  wt_val <- g_lambda_val - g_val
  return(wt_val)
}
 wts_is_mala <- numeric(length = 1e5)
 wts_is_hmc <- numeric(length = 1e5)
 
 for (i in 1:1e5) {
   wts_is_mala[i] <- wts_func(output_single_run[[1]][i,], lamb_coeff, 
                              alpha_hat, sigma2_hat, k = k, grid = x)
   wts_is_hmc[i] <- wts_func(output_single_run[[3]][i,], lamb_coeff, 
                              alpha_hat, sigma2_hat, k = k, grid = x)
 }
 
 #### Quantile for pi
 # MALA
 augm_mat <- cbind(output_single_run[[1]],wts_is_mala)
 
 upper_quant_mala_pi <- numeric(length = length(y))
 lower_quant_mala_pi <- numeric(length = length(y))
 post_med_mala_pi <- numeric(length = length(y))
 signif_level <- 0.025
 
 for (i in 1:length(y)) 
 {
   initial_mat <- quant(i, augm_mat)
   mat_sum <- apply(initial_mat, 2, sum)
   wts_prop <- initial_mat[,2]/mat_sum[2]
   final_mat <- cbind(initial_mat[1,], cumsum(wts_prop))
   lower_index <- min(which(final_mat[,2] >= signif_level))
   upper_index <- min(which(final_mat[,2] >= (1 - signif_level)))
   med_index <- min(which(final_mat[,2] >= 0.5))
   upper_quant_mala_pi[i] <- initial_mat[upper_index,1]
   lower_quant_mala_pi[i] <- initial_mat[lower_index,1]
   post_med_mala_pi[i] <- initial_mat[med_index,1]
 }
 
 ## HMC
 augm_mat <- cbind(output_single_run[[3]],wts_is_hmc)
 
 upper_quant_hmc_pi <- numeric(length = length(y))
 lower_quant_hmc_pi <- numeric(length = length(y))
 post_med_hmc_pi <- numeric(length = length(y))
 signif_level <- 0.025
 
 for (i in 1:length(y)) 
 {
   initial_mat <- quant(i, augm_mat)
   mat_sum <- apply(initial_mat, 2, sum)
   wts_prop <- initial_mat[,2]/mat_sum[2]
   final_mat <- cbind(initial_mat[1,], cumsum(wts_prop))
   lower_index <- min(which(final_mat[,2] >= signif_level))
   upper_index <- min(which(final_mat[,2] >= (1 - signif_level)))
   med_index <- min(which(final_mat[,2] >= 0.5))
   upper_quant_hmc_pi[i] <- initial_mat[upper_index,1]
   lower_quant_hmc_pi[i] <- initial_mat[lower_index,1]
   post_med_hmc_pi[i] <- initial_mat[med_index,1]
 }
 
 ##### Quantile for pi^lambda
 ## MALA
 upper_quant_pi_lambda <- numeric(length = length(y))
 lower_quant_pi_lambda <- numeric(length = length(y))
 post_med_pi_lambda <- numeric(length = length(y))
 for(i in 1:length(y))
 {
   upper_quant_pi_lambda[i] <- quantile(output_single_run[[1]][,i], probs = 0.975)
   lower_quant_pi_lambda[i] <- quantile(output_single_run[[1]][,i], probs = 0.025)
   post_med_pi_lambda[i] <- quantile(output_single_run[[1]][,i], probs = 0.5)
 }
 
 ## HMC
 upper_quant_pi_lambda.hmc <- numeric(length = length(y))
 lower_quant_pi_lambda.hmc <- numeric(length = length(y))
 post_med_pi_lambda.hmc <- numeric(length = length(y))
 for(i in 1:length(y))
 {
   upper_quant_pi_lambda.hmc[i] <- quantile(output_single_run[[3]][,i], probs = 0.975)
   lower_quant_pi_lambda.hmc[i] <- quantile(output_single_run[[3]][,i], probs = 0.025)
   post_med_pi_lambda.hmc[i] <- quantile(output_single_run[[3]][,i], probs = 0.5)
 }
 
 
 pdf(file = "tf_quantiles_mala.pdf", width = 12, height = 6)
 par(mfrow = c(1,2))

  ## pi
 post_med <- post_med_mala_pi
 upper_quant <- upper_quant_mala_pi
 lower_quant <- lower_quant_mala_pi
 dataset <- data.frame(x, y, lower_quant, upper_quant, post_med)
 plot <- ggplot(dataset, aes(x, y,group = )) + geom_point() +
   geom_line(aes(x=c(1:100), y=post_med), col = "red")
 conf_bands <- plot + geom_ribbon(aes(ymin = lower_quant, ymax = upper_quant), alpha = 0.3) +
   ggtitle("Piecewise linear model") + labs(x = "variable") + labs(y = "data")
 conf_bands

  ## pi lambda
 post_med <- post_med_pi_lambda
 upper_quant <- upper_quant_pi_lambda
 lower_quant <- lower_quant_pi_lambda
 dataset <- data.frame(x, y, lower_quant, upper_quant, post_med)
 plot <- ggplot(dataset, aes(x, y,group = )) + geom_point() +
   geom_line(aes(x=c(1:100), y=post_med), col = "red")
 conf_bands <- plot + geom_ribbon(aes(ymin = lower_quant, ymax = upper_quant), alpha = 0.3) +
   ggtitle("Piecewise linear model (pi^lambda)") + labs(x = "variable") + labs(y = "data")
 conf_bands
 dev.off()
 
 
 
 pdf(file = "tf_quantiles_hmc.pdf", width = 12, height = 6)
 par(mfrow = c(1,2))

  ## pi
 post_med <- post_med_hmc_pi
 upper_quant <- upper_quant_hmc_pi
 lower_quant <- lower_quant_hmc_pi
 dataset <- data.frame(x, y, lower_quant, upper_quant, post_med)
 plot <- ggplot(dataset, aes(x, y,group = )) + geom_point() +
   geom_line(aes(x=c(1:100), y=post_med), col = "red")
 conf_bands <- plot + geom_ribbon(aes(ymin = lower_quant, ymax = upper_quant), alpha = 0.3) +
   ggtitle("Piecewise linear model(pi)") + labs(x = "variable") + labs(y = "data")
 conf_bands

  ## pi lambda
 post_med <- post_med_pi_lambda.hmc
 upper_quant <- upper_quant_pi_lambda.hmc
 lower_quant <- lower_quant_pi_lambda.hmc
 dataset <- data.frame(x, y, lower_quant, upper_quant, post_med)
 plot <- ggplot(dataset, aes(x, y,group = )) + geom_point() +
   geom_line(aes(x=c(1:100), y=post_med), col = "red")
 conf_bands <- plot + geom_ribbon(aes(ymin = lower_quant, ymax = upper_quant), alpha = 0.3) +
   ggtitle("Piecewise linear model (pi^lambda)") + labs(x = "variable") + labs(y = "data")
 conf_bands
 dev.off()
 
 ################    ACF plots
 
pdf("acf_tf.pdf", height = 6, width = 12)
par(mfrow = c(1,2))

acf_is <- acf(output_single_run[[1]][,1], plot = FALSE)$acf
acf_pxm <- acf(output_single_run[[2]][,1], plot = FALSE)$acf

plot(1:length(acf_is), acf_is, col = "blue", type = "l",
     xlab = "Lag", ylab = "Autocorrelation", ylim = c(-0.2, 1), main = "MALA")
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
     xlab = "Lag", ylab = "Autocorrelation", ylim = c(-0.2, 1), main = "HMC")
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

