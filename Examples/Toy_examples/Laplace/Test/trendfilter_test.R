
source("TF_functions.R")
load("warmup.Rdata")

iter_mala <- 1e4
iter_hmc <- 1e3
lamb_coeff <- 0.001
D_mat <- getD(k=1, n=1e2, x)   #  D matrix
delta_samp_is <- 0.0015
delta_samp_pxm <- 0.0008


time_ismala <- system.time(mala.is <- mymala(y, alpha_hat, sigma2_hat, k=1, grid=x, iter = iter_mala, 
                                             delta = delta_samp_is, start = warmup_end_iter))

time_pxmala <- system.time(pxmala.run <- px.mala(y, alpha_hat, sigma2_hat, k=1, grid=x, iter = iter_mala, 
                                                 delta = delta_samp_pxm, start = warmup_end_iter))

# time_ishmc <- system.time(my.hmc <- myhmc(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter_hmc,
#                                           eps_hmc = 0.015, L = 100, start = warmup_end_iter))
# 
# time_pxhmc <- system.time(px.hmc <- pxhmc(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter_hmc,
#                                           eps_hmc = 0.00035, L = 100, start = warmup_end_iter))

pdf(file = "Test/sample_compar.pdf", width = 10, height = 7)
# Colours with good visibility when overlaid
col_mala_is <- adjustcolor("red", alpha.f = 0.45)
col_mala_px   <- adjustcolor("blue", alpha.f = 0.45)

# Initialise plot using both matrices to determine common y-axis limits
matplot(
  t(mat_mala_is),
  type = "l",
  lty = 1,
  col = col_mala_is,
  xlab = "time",
  ylab = "value",
  ylim = range(c(mat_mala_is, mat_mala_px)),
  main = "MALA-IS vs PxMALA"
)

# Overlay PxMALA trajectories
matlines(
  t(mat_mala_px),
  type = "l",
  lty = 1,
  col = col_mala_px
)

legend(
  "topright",
  legend = c("MALA-IS", "PxMALA"),
  col = c("red", "blue"),
  lty = 1,
  lwd = 2,
  bty = "n"
)
dev.off()
