
#####################################################################
########################## Laplace output ###########################
#####################################################################
set.seed(100)
source("laplace_functions.R")
library(mcmcse)

d <- 20
beta <- rep(1, d)
start <- rep(0.5, d)
lamb <- 1
delta_px <- .45
delta_is <- 1.2
samp <- 1e5   ##### sample size

################ functions for mixture sampling ################

samp_mix_gaussian <- function(size, means, sigma2)
{
  dim <- ncol(means)
  mean_size <- nrow(means)
  index <- sample.int(mean_size, size, replace = TRUE)
  norm_mat <- matrix(rnorm(size * dim, sd = sqrt(sigma2)),
              nrow = size, ncol = dim)
  samp_mat <- norm_mat + means[index, ]
  return(samp_mat)
}

mixture_output <- function(means, sigma2, sample_size){    ##### means denote mixture mean vector
  mc_sample <- samp_mix_gaussian(sample_size, means, sigma2) # sampling
  mean_size <- nrow(means)
  log_dens <- matrix(0, nrow = sample_size, mean_size)
  log_mix_dens <- numeric(sample_size)
  for(j in 1:mean_size){
    log_dens[,j] <- apply(mc_sample, 1, function(z) {
      sum(dnorm(z, mean = means[j,], sd = sqrt(sigma2), log = TRUE))})
  }
  for (i in 1:sample_size) {
    m <- max(log_dens[i, ])
    log_mix_dens[i] <- m + log(mean(exp(log_dens[i, ] - m)))
  }
  log_wt_pi <- apply(mc_sample, 1, function(x) {log_pi(beta, x)})
  wts_mixture <- log_wt_pi - log_mix_dens
  output <- list(mc_sample, wts_mixture)
  return(output)
}

################### For thinning ###################

iter_thin <- 1e3
pxm_for_thin <- pxmala(beta, start, lambda = lamb, iter_thin, delta_px)
mixture_means <- pxm_for_thin[seq(1, nrow(pxm_for_thin), by = 10),]  ### thinned MCMC sample
#est_covmat <- cov(pxm_for_thin). ###### For gaussian importance sample

################. For checking an optimal sigma2 ##################
start_time <- Sys.time()
iter_sigma2 <- 1e4
sigma2_vec <- seq(0.1, 3,  length = 50)
imp_ess <- numeric(length = length(sigma2_vec))
for(i in 1:length(sigma2_vec))
{
  mix_sigma <- mixture_output(mixture_means, sigma2_vec[i], sample_size = iter_sigma2)
  exp_wts <- exp(mix_sigma[[2]])
  imp_ess[i] <- (mean(exp_wts))^2 / mean(exp_wts^2)
}
end_time <- Sys.time()
plot(sigma2_vec, imp_ess, type = "l")

################ Importance proposal Gaussian mixture and MCMC. ################
 sigma2 <- 2
 mix.sample <- system.time(mixture.output <- mixture_output(mixture_means, sigma2 = sigma2,
                                                                       sample_size = samp))
 mym.sample <- system.time(mym.output <- mymala(beta, start, lambda = lamb,
                                                iter = samp, delta_is))
 
 pdf(file = "plots/weight_comp_mvl.pdf", width = 10, height = 6)
 breaks <- pretty(range(mym.output[[2]], mixture.output[[2]]), n = 300)
 h1 <- hist(mym.output[[2]], breaks = breaks, plot = FALSE)
 h2 <- hist(mixture.output[[2]], breaks = breaks, plot = FALSE)

 ylim <- c(0, max(h1$density, h2$density))

 hist(mym.output[[2]], breaks = breaks, probability = TRUE,
      xlim = range(mym.output[[2]], mixture.output[[2]]), ylim = ylim,
      col = rgb(1, 0, 0, 0.4), border = "red",
      main = "Comparison of log of importance weights",
      xlab = "Weight", ylab = "Frequency")

 hist(mixture.output[[2]], breaks = breaks, probability = TRUE,
      col = rgb(0, 0, 1, 0.4), border = "blue",
      add = TRUE)

 legend("topright",
        legend = c("ISMALA", "Gaussian Mixture"),
        fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)))
 dev.off()
