
#####################################################################
################## Laplace output #######################
#####################################################################
set.seed(100)
source("laplace_functions.R")
library(mcmcse)

d <- 1
beta <- rep(1, d)
start <- rep(0.5, d)
lamb <- 5
delta_px <- 1.85
delta_is <- 13
iter_samp <- 1e6
################### For thinning ###################

iter_thin <- 1e4
pxm_for_thin <- pxmala(beta, start, lambda = lamb, iter_thin, delta_px)
mixture_means <- pxm_for_thin[seq(1, nrow(pxm_for_thin), by = 100)]  ### thinned MCMC sample

################ Importance proposal Gaussian mixture ################
 sigma2 <- 0.01
 mc_sample <- numeric(length = iter_samp)
 mixture_index <- length(mixture_means)
 samp_index <- sample(seq(1, mixture_index), iter_samp, replace = TRUE)
 
 mixture_output <- function(means){
   mc_sample <- rnorm(n = iter_samp, mean = means[samp_index],sd = sqrt(sigma2))
   ############### weights calculation ################
   mixture_density <- vapply(mc_sample, function(z) {
     mean(dnorm(z, mean = means, sd = sqrt(sigma2)))
   }, numeric(1))
   log_wt_pi <- vapply(mc_sample, function(x) {log_pi(beta, x)}, numeric(1))
   wts_mixture <- log_wt_pi - log(mixture_density)
   output <- list(mc_sample, wts_mixture)
   return(output)
 }
 
############# For MCMC #############

 mym.sample <- system.time(mym.output <- mymala(beta, start, lambda = lamb, 
                                               iter = iter_samp, delta_is))
 mix.sample <- system.time(mixture.output <- mixture_output(mixture_means))
 
 pdf(file = "plots/weight_comp.pdf", width = 10, height = 6)
 breaks <- pretty(range(mym.output[[2]], mixture.output[[2]]), n = 300)
 h1 <- hist(mym.output[[2]], breaks = breaks, plot = FALSE)
 h2 <- hist(mixture.output[[2]], breaks = breaks, plot = FALSE)
 
 ylim <- c(0, max(h1$density, h2$density))
 
 hist(mym.output[[2]], breaks = breaks, probability = TRUE,
      xlim = range(mym.output[[2]], mixture.output[[2]]), ylim = ylim,
      col = rgb(1, 0, 0, 0.4), border = "red",
      main = "Comparison of Importance Weights",
      xlab = "Weight")
 
 hist(mixture.output[[2]], breaks = breaks, probability = TRUE,
      col = rgb(0, 0, 1, 0.4), border = "blue",
      add = TRUE)
 
 legend("topright",
        legend = c("ISMALA", "Gaussian Mixture"),
        fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)))
 dev.off()
 