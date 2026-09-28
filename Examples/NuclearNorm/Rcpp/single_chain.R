###############################
# Single chains using Rcpp
# Run from the project root:
# source("Rcpp/single_chain.R")
###############################

source("Rcpp/nuclear_norm_functions.R")
load("warmup_chain.Rdata")

# Use the original data seed, then continue its RNG stream for sampling.
data <- generate_checkerboard(seed = 8024248)
n <- data$n
checker <- data$checker
noise <- data$noise
image_mat <- data$image_mat
x <- data$x
y <- data$y

iter <- 1e4
lamb_coeff <- 1e-4
sigma2_hat <- 0.01
alpha_hat <- 1.15 / sigma2_hat
step_ismala <- 0.00012
step_pxmala <- 0.0001
step_isb <- 0.0001
step_pxb <- 0.0001
eps_is <- 0.008
eps_px <- 0.0038
L <- 10
include_barker <- FALSE
verbose <- TRUE # Ten progress updates and the final acceptance rate per sampler.

######################### MALA #########################

result_pxm <- px.mala(y = y, alpha = alpha_hat, lambda = lamb_coeff,
                     sigma2 = sigma2_hat, iter = iter,
                     delta = step_pxmala, start = warmup_end_iter, verbose = verbose)

result_ism <- mymala(y = y, alpha = alpha_hat, lambda = lamb_coeff,
                    sigma2 = sigma2_hat, iter = iter,
                    delta = step_ismala, start = warmup_end_iter, verbose = verbose)

output_single_mala <- list(result_ism[[1]], result_pxm)
save(output_single_mala, file = "Rcpp/output_single_chain_mala.Rdata")

######################### Barker #########################

if (include_barker) {
  result_pxb <- px.barker(y = y, alpha = alpha_hat, lambda = lamb_coeff,
                         sigma2 = sigma2_hat, iter = iter,
                         delta = step_pxb, start = warmup_end_iter, verbose = verbose)

  result_isb <- mybarker(y = y, alpha = alpha_hat, lambda = lamb_coeff,
                        sigma2 = sigma2_hat, iter = iter,
                        delta = step_isb, start = warmup_end_iter, verbose = verbose)

  output_single_bark <- list(result_isb[[1]], result_pxb[[1]])
  save(output_single_bark, file = "Rcpp/output_single_chain_bark.Rdata")
}

######################### HMC #########################

result_pxhmc <- pxhmc(y = y, alpha = alpha_hat, lambda = lamb_coeff,
                     sigma2 = sigma2_hat, iter = iter,
                     eps_hmc = eps_px, L = L, start = warmup_end_iter, verbose = verbose)

result_ishmc <- myhmc(y = y, alpha = alpha_hat, lambda = lamb_coeff,
                     sigma2 = sigma2_hat, iter = iter,
                     eps_hmc = eps_is, L = L, start = warmup_end_iter, verbose = verbose)

output_single_hmc <- list(result_ishmc[[1]], result_pxhmc[[1]])
save(output_single_hmc, file = "Rcpp/output_single_chain_hmc.Rdata")
