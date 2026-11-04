set.seed(100)
###########################################
# Comparing MY-IS with Mixture-based IS
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


# Parameters

d <- 20
beta <- rep(1, d)
start <- rep(0.5, d)
df <- 10

pilot_n <- 1e4
pilot_thin <- 100

sigma2_vec <- seq(0.05, 3, length.out = 50)
iter_sigma2 <- 1e4
reps <- 100
samp <- 1e4

out_dir <- "mixture_results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# 1. Laplace / Gaussian

laplace_gaussian_target <- make_laplace_target(beta)
laplace_gaussian_lambda <- 1
laplace_gaussian_delta_px <- 0.45
laplace_gaussian_delta_my <- 1.2

laplace_gaussian_pilot <- proximal_mala(
  laplace_gaussian_target, start,
  laplace_gaussian_lambda, laplace_gaussian_delta_px,
  pilot_n, mode = "px"
)

laplace_gaussian_means <- laplace_gaussian_pilot$samples[
  seq(1, pilot_n, by = pilot_thin), , drop = FALSE
]

laplace_gaussian_ess <- matrix(0, nrow = reps, ncol = length(sigma2_vec))

for (r in seq_len(reps)) 
{
  print(r)
  for (i in seq_along(sigma2_vec)) {
    laplace_gaussian_tune <- mixture_output(
      laplace_gaussian_target, laplace_gaussian_means, sigma2_vec[i],
      iter_sigma2, "gaussian", df
    )
    
    laplace_gaussian_ess[r, i] <- laplace_gaussian_tune$diagnostics$relative_ess
  }
}

laplace_gaussian_mean_ess <- colMeans(laplace_gaussian_ess)
laplace_gaussian_sigma2 <- sigma2_vec[which.max(laplace_gaussian_mean_ess)]

pdf(file.path(out_dir, "laplace_gaussian_ess.pdf"), width = 10, height = 6)

plot(
  sigma2_vec, laplace_gaussian_mean_ess,
  type = "b", pch = 19, col = "steelblue",
  xlab = expression(sigma^2), ylab = "Mean ESS / n"
)
abline(v = laplace_gaussian_sigma2, col = "firebrick", lty = 2)

dev.off()

laplace_gaussian_output <- mixture_output(
  laplace_gaussian_target, laplace_gaussian_means, laplace_gaussian_sigma2,
  samp, "gaussian", df
)

laplace_gaussian_my <- proximal_mala(
  laplace_gaussian_target, start,
  laplace_gaussian_lambda, laplace_gaussian_delta_my,
  samp, mode = "my"
)
laplace_gaussian_my$diagnostics <- weight_summary(laplace_gaussian_my$log_weights)

laplace_gaussian_mix_logs <- laplace_gaussian_output$log_weights
laplace_gaussian_my_logs <- laplace_gaussian_my$log_weights

laplace_gaussian_limits <- range(
  laplace_gaussian_mix_logs, laplace_gaussian_my_logs
)
laplace_gaussian_breaks <- pretty(laplace_gaussian_limits, n = 300)

laplace_gaussian_mix_hist <- hist(
  laplace_gaussian_mix_logs,
  breaks = laplace_gaussian_breaks, plot = FALSE
)
laplace_gaussian_my_hist <- hist(
  laplace_gaussian_my_logs,
  breaks = laplace_gaussian_breaks, plot = FALSE
)

laplace_gaussian_ylim <- c(
  0, max(laplace_gaussian_mix_hist$density, laplace_gaussian_my_hist$density)
)

pdf(file.path(out_dir, "laplace_gaussian_weights.pdf"), width = 10, height = 6)

hist(
  laplace_gaussian_my_logs, breaks = laplace_gaussian_breaks, probability = TRUE,
  xlim = laplace_gaussian_limits, ylim = laplace_gaussian_ylim,
  col = rgb(1, 0, 0, 0.4), border = "red",
  xlab = "Log importance weight", ylab = "Density"
)
hist(
  laplace_gaussian_mix_logs, breaks = laplace_gaussian_breaks, probability = TRUE,
  col = rgb(0, 0, 1, 0.4), border = "blue", add = TRUE
)
legend(
  "topright", c("MY-MALA IS", "Mixture IS"), bty = "n",
  fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4))
)

dev.off()

laplace_gaussian_summary <- data.frame(
  sigma2 = laplace_gaussian_sigma2,
  mixture_ess = laplace_gaussian_output$diagnostics$ess,
  mixture_relative_ess = laplace_gaussian_output$diagnostics$relative_ess,
  my_weight_ess = laplace_gaussian_my$diagnostics$ess,
  my_weight_relative_ess = laplace_gaussian_my$diagnostics$relative_ess
)


# 2. Laplace / Student t (df = 5)

laplace_student_target <- make_laplace_target(beta)
laplace_student_lambda <- 1
laplace_student_delta_px <- 0.45
laplace_student_delta_my <- 1.2

laplace_student_pilot <- proximal_mala(
  laplace_student_target, start,
  laplace_student_lambda, laplace_student_delta_px,
  pilot_n, mode = "px"
)

laplace_student_means <- laplace_student_pilot$samples[
  seq(1, pilot_n, by = pilot_thin), , drop = FALSE
]

laplace_student_ess <- matrix(0, nrow = reps, ncol = length(sigma2_vec))

for (r in seq_len(reps)) 
{
  print(r)
  for (i in seq_along(sigma2_vec)) {
    laplace_student_tune <- mixture_output(
      laplace_student_target, laplace_student_means, sigma2_vec[i],
      iter_sigma2, "student", df
    )
    
    laplace_student_ess[r, i] <- laplace_student_tune$diagnostics$relative_ess
  }
}

laplace_student_mean_ess <- colMeans(laplace_student_ess)
laplace_student_sigma2 <- sigma2_vec[which.max(laplace_student_mean_ess)]

pdf(file.path(out_dir, "laplace_student_ess.pdf"), width = 10, height = 6)

plot(
  sigma2_vec, laplace_student_mean_ess,
  type = "b", pch = 19, col = "steelblue",
  xlab = expression(sigma^2), ylab = "Mean ESS / n"
)
abline(v = laplace_student_sigma2, col = "firebrick", lty = 2)
dev.off()

laplace_student_output <- mixture_output(
  laplace_student_target, laplace_student_means, laplace_student_sigma2,
  samp, "student", df
)


laplace_student_my <- proximal_mala(
  laplace_student_target, start,
  laplace_student_lambda, laplace_student_delta_my,
  samp, mode = "my"
)
laplace_student_my$diagnostics <- weight_summary(laplace_student_my$log_weights)

laplace_student_mix_logs <- laplace_student_output$log_weights
laplace_student_my_logs <- laplace_student_my$log_weights

laplace_student_limits <- range(
  laplace_student_mix_logs, laplace_student_my_logs
)
laplace_student_breaks <- pretty(laplace_student_limits, n = 300)

laplace_student_mix_hist <- hist(
  laplace_student_mix_logs,
  breaks = laplace_student_breaks, plot = FALSE
)
laplace_student_my_hist <- hist(
  laplace_student_my_logs,
  breaks = laplace_student_breaks, plot = FALSE
)

laplace_student_ylim <- c(
  0, max(laplace_student_mix_hist$density, laplace_student_my_hist$density)
)

pdf(file.path(out_dir, "laplace_student_weights.pdf"), width = 10, height = 6)

hist(
  laplace_student_my_logs, breaks = laplace_student_breaks, probability = TRUE,
  xlim = laplace_student_limits, ylim = laplace_student_ylim,
  col = rgb(1, 0, 0, 0.4), border = "red",
  xlab = "Log importance weight", ylab = "Density"
)
hist(
  laplace_student_mix_logs, breaks = laplace_student_breaks, probability = TRUE,
  col = rgb(0, 0, 1, 0.4), border = "blue", add = TRUE
)
legend(
  "topright", c("MY-MALA IS", "Mixture IS"), bty = "n",
  fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4))
)

dev.off()

laplace_student_summary <- data.frame(
  sigma2 = laplace_student_sigma2,
  mixture_ess = laplace_student_output$diagnostics$ess,
  mixture_relative_ess = laplace_student_output$diagnostics$relative_ess,
  my_weight_ess = laplace_student_my$diagnostics$ess,
  my_weight_relative_ess = laplace_student_my$diagnostics$relative_ess
)


# 3. Quartic / Gaussian

exp4_gaussian_target <- make_exp4_target(d)
exp4_gaussian_lambda <- 0.1
exp4_gaussian_delta_px <- 0.09
exp4_gaussian_delta_my <- 0.30

exp4_gaussian_pilot <- proximal_mala(
  exp4_gaussian_target, start,
  exp4_gaussian_lambda, exp4_gaussian_delta_px,
  pilot_n, mode = "px"
)

exp4_gaussian_means <- exp4_gaussian_pilot$samples[
  seq(1, pilot_n, by = pilot_thin), , drop = FALSE
]

exp4_gaussian_ess <- matrix(0, nrow = reps, ncol = length(sigma2_vec))

for (r in seq_len(reps)) 
{
  print(r)
  for (i in seq_along(sigma2_vec)) {
    exp4_gaussian_tune <- mixture_output(
      exp4_gaussian_target, exp4_gaussian_means, sigma2_vec[i],
      iter_sigma2, "gaussian", df
    )
    
    exp4_gaussian_ess[r, i] <- exp4_gaussian_tune$diagnostics$relative_ess
  }
}

exp4_gaussian_mean_ess <- colMeans(exp4_gaussian_ess)
exp4_gaussian_sigma2 <- sigma2_vec[which.max(exp4_gaussian_mean_ess)]

pdf(file.path(out_dir, "exp4_gaussian_ess.pdf"), width = 10, height = 6)

plot(
  sigma2_vec, exp4_gaussian_mean_ess,
  type = "b", pch = 19, col = "steelblue",
  xlab = expression(sigma^2), ylab = "Mean ESS / n"
)
abline(v = exp4_gaussian_sigma2, col = "firebrick", lty = 2)

dev.off()


exp4_gaussian_output <- mixture_output(
  exp4_gaussian_target, exp4_gaussian_means, exp4_gaussian_sigma2,
  samp, "gaussian", df
)


exp4_gaussian_my <- proximal_mala(
  exp4_gaussian_target, start,
  exp4_gaussian_lambda, exp4_gaussian_delta_my,
  samp, mode = "my"
)
exp4_gaussian_my$diagnostics <- weight_summary(exp4_gaussian_my$log_weights)

exp4_gaussian_mix_logs <- exp4_gaussian_output$log_weights
exp4_gaussian_my_logs <- exp4_gaussian_my$log_weights

exp4_gaussian_limits <- range(
  exp4_gaussian_mix_logs, exp4_gaussian_my_logs
)
exp4_gaussian_breaks <- pretty(exp4_gaussian_limits, n = 300)

exp4_gaussian_mix_hist <- hist(
  exp4_gaussian_mix_logs,
  breaks = exp4_gaussian_breaks, plot = FALSE
)
exp4_gaussian_my_hist <- hist(
  exp4_gaussian_my_logs,
  breaks = exp4_gaussian_breaks, plot = FALSE
)

exp4_gaussian_ylim <- c(
  0, max(exp4_gaussian_mix_hist$density, exp4_gaussian_my_hist$density)
)

pdf(file.path(out_dir, "exp4_gaussian_weights.pdf"), width = 10, height = 6)

hist(
  exp4_gaussian_my_logs, breaks = exp4_gaussian_breaks, probability = TRUE,
  xlim = exp4_gaussian_limits, ylim = exp4_gaussian_ylim,
  col = rgb(1, 0, 0, 0.4), border = "red",
  xlab = "Log importance weight", ylab = "Density"
)
hist(
  exp4_gaussian_mix_logs, breaks = exp4_gaussian_breaks, probability = TRUE,
  col = rgb(0, 0, 1, 0.4), border = "blue", add = TRUE
)
legend(
  "topright", c("MY-MALA IS", "Mixture IS"), bty = "n",
  fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4))
)

dev.off()

exp4_gaussian_summary <- data.frame(
  sigma2 = exp4_gaussian_sigma2,
  mixture_ess = exp4_gaussian_output$diagnostics$ess,
  mixture_relative_ess = exp4_gaussian_output$diagnostics$relative_ess,
  my_weight_ess = exp4_gaussian_my$diagnostics$ess,
  my_weight_relative_ess = exp4_gaussian_my$diagnostics$relative_ess
)


# 4. Quartic / Student t (df = 5)

exp4_student_target <- make_exp4_target(d)
exp4_student_lambda <- 0.1
exp4_student_delta_px <- 0.09
exp4_student_delta_my <- 0.30


exp4_student_pilot <- proximal_mala(
  exp4_student_target, start,
  exp4_student_lambda, exp4_student_delta_px,
  pilot_n, mode = "px"
)

exp4_student_means <- exp4_student_pilot$samples[
  seq(1, pilot_n, by = pilot_thin), , drop = FALSE
]


exp4_student_ess <- matrix(0, nrow = reps, ncol = length(sigma2_vec))

for (r in seq_len(reps))
{
  print(r)
  for (i in seq_along(sigma2_vec)) {
    exp4_student_tune <- mixture_output(
      exp4_student_target, exp4_student_means, sigma2_vec[i],
      iter_sigma2, "student", df
    )
    
    exp4_student_ess[r, i] <- exp4_student_tune$diagnostics$relative_ess
  }
}

exp4_student_mean_ess <- colMeans(exp4_student_ess)
exp4_student_sigma2 <- sigma2_vec[which.max(exp4_student_mean_ess)]

pdf(file.path(out_dir, "exp4_student_ess.pdf"), width = 10, height = 6)

plot(
  sigma2_vec, exp4_student_mean_ess,
  type = "b", pch = 19, col = "steelblue",
  xlab = expression(sigma^2), ylab = "Mean ESS / n"
)
abline(v = exp4_student_sigma2, col = "firebrick", lty = 2)

dev.off()

exp4_student_output <- mixture_output(
  exp4_student_target, exp4_student_means, exp4_student_sigma2,
  samp, "student", df
)

exp4_student_my <- proximal_mala(
  exp4_student_target, start,
  exp4_student_lambda, exp4_student_delta_my,
  samp, mode = "my"
)
exp4_student_my$diagnostics <- weight_summary(exp4_student_my$log_weights)

exp4_student_mix_logs <- exp4_student_output$log_weights
exp4_student_my_logs <- exp4_student_my$log_weights

exp4_student_limits <- range(
  exp4_student_mix_logs, exp4_student_my_logs
)
exp4_student_breaks <- pretty(exp4_student_limits, n = 300)

exp4_student_mix_hist <- hist(
  exp4_student_mix_logs,
  breaks = exp4_student_breaks, plot = FALSE
)
exp4_student_my_hist <- hist(
  exp4_student_my_logs,
  breaks = exp4_student_breaks, plot = FALSE
)

exp4_student_ylim <- c(
  0, max(exp4_student_mix_hist$density, exp4_student_my_hist$density)
)

pdf(file.path(out_dir, "exp4_student_weights.pdf"), width = 10, height = 6)

hist(
  exp4_student_my_logs, breaks = exp4_student_breaks, probability = TRUE,
  xlim = exp4_student_limits, ylim = exp4_student_ylim,
  col = rgb(1, 0, 0, 0.4), border = "red",
  xlab = "Log importance weight", ylab = "Density"
)
hist(
  exp4_student_mix_logs, breaks = exp4_student_breaks, probability = TRUE,
  col = rgb(0, 0, 1, 0.4), border = "blue", add = TRUE
)
legend(
  "topright", c("MY-MALA IS", "Mixture IS"), bty = "n",
  fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4))
)

dev.off()

exp4_student_summary <- data.frame(
  sigma2 = exp4_student_sigma2,
  mixture_ess = exp4_student_output$diagnostics$ess,
  mixture_relative_ess = exp4_student_output$diagnostics$relative_ess,
  my_weight_ess = exp4_student_my$diagnostics$ess,
  my_weight_relative_ess = exp4_student_my$diagnostics$relative_ess
)
