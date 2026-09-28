# Run scripts from the project root (the folder containing warmup_chain.Rdata).
# Sourcing this file compiles the samplers; it does not run a simulation.
library(Rcpp)
library(RcppArmadillo)

Rcpp::sourceCpp(
  "Rcpp/nuclear_norm_functions.cpp",
  env = environment(),
  cacheDir = getOption("nuclear_norm.rcpp.cache_dir", file.path(tempdir(), "nuclear_norm_rcpp"))
)

# Standalone mathematical helpers remain in R. The sampler engine uses private
# C++ equivalents to avoid R callbacks inside the sampling loops.
matrix_dimension <- function(x) {
  n <- sqrt(length(x))
  if (!length(x) || n != floor(n) || any(!is.finite(x)))
    stop("Input must be a finite, nonempty vector of square length.")
  as.integer(n)
}
nucl_norm <- function(vect) {
  n <- matrix_dimension(vect)
  sum(svd(matrix(vect, n, n), nu = 0, nv = 0)$d)
}
log_pi <- function(x, y, sigma2, alpha)
  -alpha * nucl_norm(x) - sum((y - x)^2) / (2 * sigma2)
log_pilambda <- function(eta, x, lambda, y, sigma2, alpha)
  log_pi(eta, y, sigma2, alpha) - sum((eta - x)^2) / (2 * lambda)
softthreshold <- function(u, lambda) sign(u) * pmax(abs(u) - lambda, 0)
prox_func <- function(x, lambda, y, sigma2, alpha) {
  n <- matrix_dimension(x)
  if (length(x) != length(y) || any(!is.finite(y)) ||
      length(lambda) != 1L || !is.finite(lambda) || lambda <= 0 ||
      length(sigma2) != 1L || !is.finite(sigma2) || sigma2 <= 0 ||
      length(alpha) != 1L || !is.finite(alpha) || alpha < 0)
    stop("Check y, lambda, sigma2 and alpha: equal lengths, finite values, positive scales, alpha >= 0.")
  fit <- svd(matrix((lambda * y + sigma2 * x) / (lambda + sigma2), n, n))
  s <- softthreshold(fit$d, alpha * sigma2 * lambda / (lambda + sigma2))
  as.vector(fit$u %*% (s * t(fit$v)))
}
grad_logpiLam <- function(x, lambda, y, sigma2, alpha)
  (prox_func(x, lambda, y, sigma2, alpha) - x) / lambda
bark.prop <- function(x, alpha, lambda, y, sigma2, delta) {
  z <- sqrt(delta) * rnorm(length(x))
  prob <- plogis(z * grad_logpiLam(x, lambda, y, sigma2, alpha))
  x + ifelse(runif(length(x)) <= prob, z, -z)
}
log_bark.dens <- function(curr_point, prop_point, grad_curr_point, delta) {
  diff <- prop_point - curr_point
  z <- -grad_curr_point * diff
  sum(dnorm(diff, 0, sqrt(delta), log = TRUE)) - sum(pmax(z, 0) + log1p(exp(-abs(z))))
}

positive_integer <- function(x, name) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      x < 1 || x != floor(x) || x > .Machine$integer.max)
    stop(name, " must be a positive integer.")
  as.integer(x)
}

mymala <- function(y, alpha, lambda, sigma2, iter, delta, start, verbose = TRUE)
  nn_sample_cpp(y, alpha, lambda, sigma2, positive_integer(iter, "iter"), delta,
                1L, start, "mala", TRUE, verbose)
px.mala <- function(y, alpha, lambda, sigma2, iter, delta, start, verbose = TRUE)
  nn_sample_cpp(y, alpha, lambda, sigma2, positive_integer(iter, "iter"), delta,
                1L, start, "mala", FALSE, verbose)
mybarker <- function(y, alpha, lambda, sigma2, iter, delta, start, verbose = TRUE)
  nn_sample_cpp(y, alpha, lambda, sigma2, positive_integer(iter, "iter"), delta,
                1L, start, "barker", TRUE, verbose)
px.barker <- function(y, alpha, lambda, sigma2, iter, delta, start, verbose = TRUE)
  nn_sample_cpp(y, alpha, lambda, sigma2, positive_integer(iter, "iter"), delta,
                1L, start, "barker", FALSE, verbose)
myhmc <- function(y, alpha, lambda, sigma2, iter, eps_hmc, L, start, verbose = TRUE)
  nn_sample_cpp(y, alpha, lambda, sigma2, positive_integer(iter, "iter"), eps_hmc,
                positive_integer(L, "L"), start, "hmc", TRUE, verbose)
pxhmc <- function(y, alpha, lambda, sigma2, iter, eps_hmc, L, start, verbose = TRUE)
  nn_sample_cpp(y, alpha, lambda, sigma2, positive_integer(iter, "iter"), eps_hmc,
                positive_integer(L, "L"), start, "hmc", FALSE, verbose)

# Original checkerboard construction and data-generation seed.
generate_checkerboard <- function(n = 64L, a = 8L, seed = 8024248L) {
  n <- positive_integer(n, "n")
  a <- positive_integer(a, "a")
  if (n %% (2 * a) != 0) stop("n must be a multiple of 2*a.")
  set.seed(seed)

  vec.mat <- rep(c(1, 0), each = a, times = n/(2*a))
  checker <- matrix(0, nrow = n, ncol = n)
  for(j in seq(1, n/2, by = a))
  {
    for(k in 1:a)
    {
      checker[j+k-1, ] <- vec.mat
    }
    vec.mat <- rev(vec.mat)
  }
  for(j in seq((n/2+1), n, by = a))
  {
    for(k in 1:a)
    {
      checker[j+k-1, ] <- (vec.mat == 0)*(vec.mat) + (vec.mat == 1)*(vec.mat - .50)
    }
    vec.mat <- rev(vec.mat)
  }
  noise <- matrix(rnorm(n^2, 0, sqrt(0.01)), nrow = n, ncol = n)
  image_mat <- checker + noise

  # as.vector() uses the same column ordering as vec() in the original code.
  x <- as.vector(checker)
  y <- as.vector(image_mat)
  list(n = n, checker = checker, noise = noise, image_mat = image_mat, x = x, y = y)
}

# Retain mcmcse's batch-means covariance estimator in R. The delta method for
# sum(X*w)/sum(w) uses the tuple (X*w, w), not (X*w, exp(w)).
asymp_cov_func <- function(chain, weights) {
  if (!requireNamespace("mcmcse", quietly = TRUE)) stop("Install the mcmcse package.")
  if (!is.matrix(chain) || nrow(chain) != length(weights) ||
      any(!is.finite(chain)) || any(!is.finite(weights)) ||
      any(weights < 0) || !any(weights > 0))
    stop("Supply a finite chain and one finite, nonnegative weight per row, with positive sum.")
  weights <- weights / max(weights)
  wts_mean <- mean(weights)
  vapply(seq_len(ncol(chain)), function(j) {
    num <- chain[, j] * weights
    is_est <- sum(num) / sum(weights)
    Sigma <- mcmcse::mcse.multi(cbind(num, weights))$cov
    derivative <- c(1, -is_est) / wts_mean
    as.numeric(crossprod(derivative, Sigma %*% derivative))
  }, numeric(1))
}

importance_weights <- function(log_weights) {
  if (!length(log_weights) || any(!is.finite(log_weights)))
    stop("Log weights must be nonempty and finite.")
  exp(log_weights - max(log_weights))
}
