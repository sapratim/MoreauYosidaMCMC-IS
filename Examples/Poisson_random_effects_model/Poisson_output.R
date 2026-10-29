################################################################################
################## Poisson example output visualisation ########################
################################################################################

################ ACF plots ################
dim <- I + 1
lag.max <- 50

################## MALA  #####################
load("Rcpp/output_single_mala.Rdata")
acf_ism <- acf(output_chain_mala[[1]][,1], plot = FALSE, lag.max = lag.max)$acf
acf_pxm <- acf(output_chain_mala[[2]][,1], plot = FALSE, lag.max = lag.max)$acf

diff.acf <- matrix(0, ncol = dim, nrow = lag.max + 1)
diff.acf[,1] <- acf_ism - acf_pxm

for (i in 2:dim) 
{
  acf_ism <- acf(output_chain_mala[[1]][,i], plot = FALSE,lag.max = lag.max)$acf
  acf_pxm <- acf(output_chain_mala[[2]][,i], plot = FALSE,lag.max = lag.max)$acf
  diff.acf[,i] <- acf_ism - acf_pxm
}

# Make boxplot of ACFs
pdf("plots/poisson_acf_MALA.pdf", height = 6, width = 6)
par(mfrow = c(1,1))
boxplot(t(diff.acf),
        xlab = "Lags", col = "pink",
        ylab = "Difference in ACF of MALAs",ylim = c(-1.3, 0.5),
        names = 0:lag.max, show.names = TRUE)

dev.off()

################## True Barker  #####################
load("Rcpp/output_single_bark.Rdata")
acf_isb <- acf(output_chain_bark[[1]][,1], plot = FALSE, lag.max = lag.max)$acf
acf_trub <- acf(output_chain_bark[[3]][,1], plot = FALSE, lag.max = lag.max)$acf

diff.acf <- matrix(0, ncol = dim, nrow = lag.max + 1)
diff.acf[,1] <- acf_isb - acf_trub

for (i in 2:dim) 
{
  acf_isb <- acf(output_chain_bark[[1]][,i], plot = FALSE, lag.max = lag.max)$acf
  acf_trub <- acf(output_chain_bark[[3]][,i], plot = FALSE, lag.max = lag.max)$acf
  diff.acf[,i] <- acf_isb - acf_trub
}

# Make boxplot of ACFs
pdf("plots/poisson_acf_Bark.pdf", height = 6, width = 6)
par(mfrow = c(1,1))
boxplot(t(diff.acf),
        xlab = "Lags", col = "pink",
        ylab = "Difference in ACF of Barkers",ylim = c(-1.3, 0.5),
        names = 0:lag.max, show.names = TRUE)

dev.off()



######################## HMC ########################
load("Rcpp/output_single_hmc.Rdata")



acf_is_hmc <- acf(output_chain_hmc[[1]][,1], plot = FALSE, lag.max = lag.max)$acf
acf_pxhmc <- acf(output_chain_hmc[[2]][,1], plot = FALSE, lag.max = lag.max)$acf

diff.acf <- matrix(0, ncol = dim, nrow = lag.max + 1)
diff.acf[,1] <- acf_is_hmc - acf_pxhmc

for (i in 2:dim) 
{
  acf_is_hmc <- acf(output_chain_hmc[[1]][,i], plot = FALSE, lag.max = lag.max)$acf
  acf_pxhmc <- acf(output_chain_hmc[[2]][,i], plot = FALSE, lag.max = lag.max)$acf
  diff.acf[,i] <- acf_is_hmc - acf_pxhmc
}

pdf("plots/poisson_acf_HMC.pdf", height = 6, width = 6)
par(mfrow = c(1,1))
# Make boxplot of ACFs
boxplot(t(diff.acf),
        xlab = "Lags", col = "pink",
        ylab = "Difference in ACF of HMCs",ylim = c(-1.3, 0.5),
        names = 0:lag.max, show.names = TRUE)

dev.off()

##################  Boxplots of marginal efficiency ##################

load("Rcpp/output_poisson.Rdata")
reps <- length(output_poisson)

#### Marginal variance comparison MALA
#x <- c(seq_len(reps))
margvar_ism <- matrix(0, nrow = reps, ncol = dim)
margvar_pxm <- matrix(0, nrow = reps, ncol = dim)
for (i in seq_len(reps))
{
  asympmat_ism <- matrix(unlist(output_poisson[[i]][[1]]), nrow = dim, ncol = dim, byrow = T)
  asympmat_pxm <- matrix(unlist(output_poisson[[i]][[2]]), nrow = dim, ncol = dim, byrow = T)
  for (j in 1:dim)
  {
    margvar_ism[i,j] <- asympmat_ism[j,j]
    margvar_pxm[i,j] <- asympmat_pxm[j,j]
  }
}
time_ism <- sapply(output_poisson, function(x) x[[8]])
time_pxm <- sapply(output_poisson, function(x) x[[9]])  
  
rel_eff_mat_mala <- (margvar_pxm*time_pxm)/(margvar_ism*time_ism)

###### Component-wise average relative efficiency  ######
mu_mala <- colMeans(rel_eff_mat_mala, na.rm = TRUE)
lower_mala <- mu_mala - 1.96*apply(rel_eff_mat_mala, 2, sd)/sqrt(100)
upper_mala <- mu_mala + 1.96*apply(rel_eff_mat_mala, 2, sd)/sqrt(100)


pdf(file = "plots/poisson_unc_mala.pdf", width = 10, height = 5)
index <- 1:ncol(rel_eff_mat_mala)
plot(index, mu_mala, pch = 16, col = "pink", ylim = c(.5, 10),
     xlab = "Component", ylab = "Average efficiency", cex.lab = 1.5)
segments(index, lower_mala, index, mu_mala)
segments(index, upper_mala, index, mu_mala)
dev.off()


#### Marginal variance comparison for Barker
margvar_trub <- matrix(0, nrow = reps, ncol = dim)
for (i in seq_len(reps))
{
  asympmat_trub <- matrix(unlist(output_poisson[[i]][[5]]), nrow = dim, ncol = dim, byrow = T)
  for (j in 1:dim)
  {
    margvar_trub[i,j] <- asympmat_trub[j,j]
  }
}

time_trub <- sapply(output_poisson, function(x) x[[12]])
rel_eff_mat_trubark <- (margvar_trub*time_trub)/(margvar_ism*time_ism)

mu_bark <- colMeans(rel_eff_mat_trubark, na.rm = TRUE)
lower_bark  <- mu_bark  - 1.96*apply(rel_eff_mat_trubark, 2, sd)/sqrt(100)
upper_bark  <- mu_bark  + 1.96*apply(rel_eff_mat_trubark, 2, sd)/sqrt(100)

pdf(file = "plots/poisson_unc_bark.pdf", width = 10, height = 5)
index <- 1:ncol(rel_eff_mat_mala)
plot(index, mu_bark, pch = 16, col = "pink", ylim = c(0, 8),
     xlab = "Component", ylab = "Average efficiency", cex.lab = 1.5)
segments(index, lower_bark, index, mu_bark)
segments(index, upper_bark, index, mu_bark)
dev.off() 


#### Marginal variance comparison for HMC
margvar_ish <- matrix(0, nrow = reps, ncol = dim)
margvar_pxh <- matrix(0, nrow = reps, ncol = dim)
for (i in seq_len(reps))
{
  asympmat_ish <- matrix(unlist(output_poisson[[i]][6]), nrow = dim, ncol = dim, byrow = T)
  asympmat_pxh <- matrix(unlist(output_poisson[[i]][7]), nrow = dim, ncol = dim, byrow = T)
  for (j in 1:dim)
  {
    margvar_ish[i,j] <- asympmat_ish[j,j]
    margvar_pxh[i,j] <- asympmat_pxh[j,j]
  }
}
time_ish <- sapply(output_poisson, function(x) x[[13]])
time_pxh <- sapply(output_poisson, function(x) x[[14]])  

rel_eff_mat_hmc <- (margvar_pxh*time_pxh)/(margvar_ish*time_ish)


mu_hmc <- colMeans(rel_eff_mat_hmc, na.rm = TRUE)
lower_hmc  <- mu_hmc  - 1.96*apply(rel_eff_mat_hmc, 2, sd)/sqrt(100)
upper_hmc  <- mu_hmc  + 1.96*apply(rel_eff_mat_hmc, 2, sd)/sqrt(100)

pdf(file = "plots/poisson_unc_hmc.pdf", width = 10, height = 5)
index <- 1:ncol(rel_eff_mat_mala)
plot(index, mu_hmc, pch = 16, col = "pink", ylim = c(0, 300),
     xlab = "Component", ylab = "Average efficiency", cex.lab = 1.5)
segments(index, lower_hmc, index, mu_hmc)
segments(index, upper_hmc, index, mu_hmc)
dev.off() 


###############  Histogram  ###############

avg_rel_eff_mala <- apply(rel_eff_mat_mala, 2, mean)
avg_rel_eff_hmc <- apply(rel_eff_mat_hmc, 2, mean)
avg_rel_eff_trub <-  apply(rel_eff_mat_trubark, 2, mean)

avg_rel_effs <- cbind(avg_rel_eff_mala, avg_rel_eff_trub)
colnames(avg_rel_effs) <- c("MALAs",  "Barker vs IS-MALA")

pdf(file = "plots/poisson_eff_mala-barker.pdf", height = 5, width = 6)
boxplot(avg_rel_effs, ylab = "Relative efficiency", show.names = TRUE,
  boxwex = .5, col = "pink", horizontal  = TRUE)
dev.off()

pdf(file = "plots/poisson_eff_hmc.pdf", height = 5, width = 6)
boxplot(avg_rel_eff_hmc, col = "pink", horizontal = TRUE, boxwex = .5, show.names = TRUE, 
  names = "HMCs", ylab = "Relative efficiency")
dev.off()



############## varying lambda output ################

lambda <- seq(0.0001, 0.004, length = 10)
load("Rcpp/output_variable_lambda.Rdata")

################################ MALA ###############################

matrices_mala <- lapply(output_lambda, '[[', 1)
se_mala_mat <- sapply(matrices_mala, function(x) rowMeans(x)) #avg eff over components for all lambda
mean_val_mala <- rowMeans(se_mala_mat)
se_val_mala <- apply(se_mala_mat, 1, sd) / sqrt(ncol(se_mala_mat))

################################ HMC #################################

matrices_hmc <- lapply(output_lambda, '[[', 3)
se_hmc_mat <- sapply(matrices_hmc, function(x) rowMeans(x)) #avg eff over components for all lambda
mean_val_hmc <- rowMeans(se_hmc_mat)
se_val_hmc <- apply(se_hmc_mat, 1, sd) / sqrt(ncol(se_hmc_mat))


################################ Barker #################################
matrices_bark <- lapply(output_lambda, '[[', 2)
se_bark_mat <- sapply(matrices_bark, function(x) rowMeans(x)) #avg eff over components for all lambda
mean_val_bark <- rowMeans(se_bark_mat)
se_val_bark <- apply(se_bark_mat, 1, sd) / sqrt(ncol(se_bark_mat))

pdf("plots/varying_lambda.pdf", width = 8, height = 6)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(lambda, mean_val_mala,
     type = "l",
     col = "red", log = "y",
     xlab = expression(lambda),
     ylab = "Log average relative efficiency ",
     ylim = range(c(mean_val_mala, mean_val_bark, mean_val_hmc)),
     cex.axis = 1.5, cex.lab = 1.5)
segments(lambda, mean_val_mala - 1.96*se_val_mala,
         lambda, mean_val_mala + 1.96*se_val_mala,
         col = "red", lwd = 2)

lines(lambda, mean_val_bark,
      col = "blue", type = "l")
segments(lambda, mean_val_bark - 1.96*se_val_bark,
         lambda, mean_val_bark + 1.96*se_val_bark,
         col = "blue", lwd = 2)

lines(lambda, mean_val_hmc,
      col = "orange", type = "l")
segments(lambda, mean_val_hmc - 1.96*se_val_hmc,
         lambda, mean_val_hmc + 1.96*se_val_hmc,
         col = "orange", lwd = 2)

abline(v = 0.001, lty = 2)

legend("topright", bty = "n",
       legend = c("ISMALA vs PxMALA",
                  "ISMALA vs Barker",
                  "ISHMC vs PxHMC"),
       col = c("red", "blue", "orange"),
       lty = 1, cex = 0.72)

dev.off()

