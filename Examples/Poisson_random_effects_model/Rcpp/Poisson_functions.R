# Run from the Poisson_random_effects_model project root.
# All sampling loops, targets, gradients and proximal solves execute in C++.
library(mcmcse)
library(foreach)
library(doParallel)
Rcpp::sourceCpp("Rcpp/Poisson_functions.cpp", env = environment(),
  cacheDir = getOption("poisson.rcpp.cache_dir", file.path(tempdir(), "poisson_rcpp")))

# Preserve the original data and RNG stream exactly.
set.seed(8024248)
ni_s <- 5
I <- 50
c <- 10
sigma_eta <- 3
one_mat <- rep(1, I)
data <- matrix(0, nrow = I, ncol = ni_s)
identity_mat <- diag(1, I, I)
mu <- rnorm(1, 0, c)
eta_vec <- rnorm(I, mu, sigma_eta)
tol_nr <- 1e-8
for (j in 1:ni_s) data[, j] <- rpois(I, exp(eta_vec))

# Public mathematical helpers retain their original positional arguments.
log_p <- function(eta, mu, data, sigma = sigma_eta, prior_sd = c) {
  -sum((eta - mu)^2) / (2 * sigma^2) - mu^2 / (2 * prior_sd^2) +
    sum(rowSums(data) * eta) - ncol(data) * sum(exp(eta))
}
log_plam <- function(eta, mu, data, lambda, y, sigma = sigma_eta, prior_sd = c)
  log_p(eta, mu, data, sigma, prior_sd) - sum((y - c(eta, mu))^2) / (2 * lambda)
grad_logp <- function(mu, sigma, eta, data = get("data", envir = environment(grad_logp)), prior_sd = c) {
  c(-(eta - mu) / sigma^2 - ncol(data) * exp(eta) + rowSums(data),
    sum(eta - mu) / sigma^2 - mu / prior_sd^2)
}
true_hessian <- function(sigma, eta, data = get("data", envir = environment(true_hessian)))
  -1 / sigma^2 - ncol(data) * exp(eta)
proxfunc <- function(eta, mu, lambda, eta_initial, mu_initial, sigma,
                     data = get("data", envir = environment(proxfunc)), prior_sd = c) {
  # Initial guesses are retained for API compatibility. The deterministic C++
  # initializer avoids history-dependent approximate forces in HMC.
  poisson_prox_cpp(c(eta, mu), data, lambda, sigma, prior_sd, tol_nr)
}
grad_logplam <- function(eta, mu, lambda, eta_initial, mu_initial, sigma,
                         data = get("data", envir = environment(grad_logplam)), prior_sd = c)
  (proxfunc(eta, mu, lambda, eta_initial, mu_initial, sigma, data, prior_sd) - c(eta, mu)) / lambda
bark.prop <- function(eta, mu, lambda, eta_initial, mu_initial, sigma, delta) {
  z <- sqrt(delta) * rnorm(length(eta) + 1)
  prob <- plogis(z * grad_logplam(eta, mu, lambda, eta_initial, mu_initial, sigma))
  c(eta, mu) + ifelse(runif(length(z)) <= prob, z, -z)
}
bark.prop_true <- function(eta, mu, sigma, delta) {
  z <- sqrt(delta) * rnorm(length(eta) + 1)
  prob <- plogis(z * grad_logp(mu, sigma, eta))
  c(eta, mu) + ifelse(runif(length(z)) <= prob, z, -z)
}
log_bark.dens <- function(curr_point, prop_point, grad_curr_point, delta) {
  diff <- prop_point - curr_point
  z <- -grad_curr_point * diff
  # Omit d*log(2), as in the original; it cancels in every MH ratio.
  sum(dnorm(diff, 0, sqrt(delta), log = TRUE)) - sum(pmax(z, 0) + log1p(exp(-abs(z))))
}
positive_integer <- function(x, name) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x < 1 ||
      x != floor(x) || x > .Machine$integer.max) stop(name, " must be a positive integer.")
  as.integer(x)
}

mymala <- function(eta_start, mu_start, lambda, sigma, iter, delta, data, verbose = TRUE, prior_sd = c)
  poisson_sample_cpp(c(eta_start, mu_start), data, lambda, sigma, prior_sd,
    positive_integer(iter, "iter"), delta, 1L, "mala", TRUE, FALSE, verbose, tol_nr)
px.mala <- function(eta_start, mu_start, lambda, sigma, iter, delta, data, verbose = TRUE, prior_sd = c)
  poisson_sample_cpp(c(eta_start, mu_start), data, lambda, sigma, prior_sd,
    positive_integer(iter, "iter"), delta, 1L, "mala", FALSE, FALSE, verbose, tol_nr)
mybarker <- function(eta_start, mu_start, lambda, sigma, iter, delta, data, verbose = TRUE, prior_sd = c)
  poisson_sample_cpp(c(eta_start, mu_start), data, lambda, sigma, prior_sd,
    positive_integer(iter, "iter"), delta, 1L, "barker", TRUE, FALSE, verbose, tol_nr)
px.barker <- function(eta_start, mu_start, lambda, sigma, iter, delta, data, verbose = TRUE, prior_sd = c)
  poisson_sample_cpp(c(eta_start, mu_start), data, lambda, sigma, prior_sd,
    positive_integer(iter, "iter"), delta, 1L, "barker", FALSE, FALSE, verbose, tol_nr)
barker <- function(eta_start, mu_start, sigma, iter, delta, data, verbose = TRUE, prior_sd = c)
  poisson_sample_cpp(c(eta_start, mu_start), data, 1, sigma, prior_sd,
    positive_integer(iter, "iter"), delta, 1L, "barker", FALSE, TRUE, verbose, tol_nr)
myhmc <- function(eta_start, mu_start, lambda, sigma, iter, data, eps_hmc, L, verbose = TRUE, prior_sd = c)
  poisson_sample_cpp(c(eta_start, mu_start), data, lambda, sigma, prior_sd,
    positive_integer(iter, "iter"), eps_hmc, positive_integer(L, "L"), "hmc", TRUE, FALSE, verbose, tol_nr)
pxhmc <- function(eta_start, mu_start, lambda, sigma, iter, data, eps_hmc, L, verbose = TRUE, prior_sd = c)
  poisson_sample_cpp(c(eta_start, mu_start), data, lambda, sigma, prior_sd,
    positive_integer(iter, "iter"), eps_hmc, positive_integer(L, "L"), "hmc", FALSE, FALSE, verbose, tol_nr)

importance_weights <- function(log_weights) {
  if (!length(log_weights) || anyNA(log_weights) || any(log_weights == Inf) ||
      !any(is.finite(log_weights))) stop("Log weights must be finite or -Inf, with at least one finite value.")
  exp(log_weights - max(log_weights))
}
asymp_covmat_fn <- function(chain, weights) {
  if (!is.matrix(chain) || nrow(chain) != length(weights) || any(!is.finite(chain)) ||
      any(!is.finite(weights)) || any(weights < 0) || !any(weights > 0))
    stop("Supply a finite chain and one nonnegative finite weight per row, with positive sum.")
  weights <- weights / max(weights)
  wts_mean <- mean(weights)
  num <- chain * weights
  is_est <- colSums(num) / sum(weights)
  Sigma_mat <- mcmcse::mcse.multi(cbind(num, weights))$cov
  derivative <- cbind(diag(1 / wts_mean, ncol(chain)), -is_est / wts_mean)
  derivative %*% Sigma_mat %*% t(derivative)
}
