###########################################
# mixture_functions.R
# Functions only (original code, unchanged).
# Nothing is run and nothing is saved here.
###########################################

# The two target distributions
make_laplace_target <- function(beta) {
  list(
    name = "laplace",
    d = length(beta),
    log_pi = function(x) -sum(beta * abs(x)),
    log_pi_rows = function(X) -as.vector(abs(X) %*% beta),
    prox = function(x, lambda) {
      sign(x) * pmax(abs(x) - lambda * beta, 0)
    }
  )
}

make_exp4_target <- function(d) {
  list(
    name = "exp4",
    d = d,
    log_pi = function(x) -sum(x^4),
    log_pi_rows = function(X) -rowSums(X^4),
    prox = function(x, lambda) {
      a <- sqrt(3 * lambda)
      sinh(asinh(3 * a * x) / 3) / a
    }
  )
}

# A generic proximal mala function
proximal_mala <- function(target, start, lambda, delta, n, mode) {
  state <- function(x) {
    p <- target$prox(x, lambda)
    log_pi <- target$log_pi(x)
    log_pi_lambda <- target$log_pi(p) - sum((x - p)^2) / (2 * lambda)
    
    list(
      x = x,
      drift = (p - x) / lambda,
      log_pi = log_pi,
      log_density = if (mode == "px") log_pi else log_pi_lambda
    )
  }
  
  current <- state(start)
  samples <- matrix(0, nrow = n, ncol = target$d)
  log_weights <- numeric(n)
  accepted <- 0
  
  samples[1, ] <- current$x
  log_weights[1] <- current$log_pi - current$log_density
  
  for (i in 2:n) {
    forward_mean <- current$x + (delta / 2) * current$drift
    proposal <- state(forward_mean + sqrt(delta) * rnorm(target$d))
    reverse_mean <- proposal$x + (delta / 2) * proposal$drift
    
    log_ratio <- proposal$log_density - current$log_density +
      (sum((proposal$x - forward_mean)^2) -
         sum((current$x - reverse_mean)^2)) / (2 * delta)
    
    if (log(runif(1)) <= log_ratio) {
      current <- proposal
      accepted <- accepted + 1
    }
    
    samples[i, ] <- current$x
    log_weights[i] <- current$log_pi - current$log_density
  }
  
  list(
    samples = samples,
    log_weights = log_weights,
    acceptance = accepted / (n - 1)
  )
}

# sampler that samples from a mixture distribution
# iid samples -- either Gaussian or t
sample_mixture <- function(n, means, sigma2, family, df = 5) 
{
  index <- sample.int(nrow(means), n, replace = TRUE)
  noise <- matrix(rnorm(n * ncol(means)), nrow = n)
  
  if (family == "student") {
    noise <- noise / sqrt(rchisq(n, df) / df)
  }
  means[index, , drop = FALSE] + sqrt(sigma2) * noise
}

# helper function
row_log_sum_exp <- function(X) {
  row_max <- X[cbind(seq_len(nrow(X)), max.col(X, ties.method = "first"))]
  row_max + log(rowSums(exp(X - row_max)))
}

# evaluate the density of a the log_mixture distribution
log_mixture <- function(X, means, sigma2, family, df = 5)
{
  d <- ncol(X)
  
  distance2 <- outer(rowSums(X^2), rowSums(means^2), "+") -
    2 * tcrossprod(X, means)
  
  distance2 <- pmax(distance2, 0)
  
  if (family == "gaussian")
  {
    component_logs <- -d / 2 * log(2 * pi * sigma2) -
      distance2 / (2 * sigma2)
  } else
  {
    component_logs <- lgamma((df + d) / 2) - lgamma(df / 2) -
      d / 2 * log(df * pi * sigma2) -
      (df + d) / 2 * log1p(distance2 / (df * sigma2))
  }
  
  row_log_sum_exp(component_logs) - log(nrow(means))
}

# function summarizes weights
weight_summary <- function(log_weights) 
{
  shifted <- log_weights - max(log_weights)
  weights <- exp(shifted)
  ess <- sum(weights)^2 / sum(weights^2)
  
  list(
    ess = ess,
    relative_ess = ess / length(weights),
    weights = weights / sum(weights)
  )
}


mixture_output <- function(target, means, sigma2, n, family, df = 5) 
{
  samples <- sample_mixture(n, means, sigma2, family, df)
  log_weights <- target$log_pi_rows(samples) -
    log_mixture(samples, means, sigma2, family, df)
  
  list(
    samples = samples,
    log_weights = log_weights,
    diagnostics = weight_summary(log_weights)
  )
}
