set.seed(100)
###########################################
# mixture_run.R
# Comparing MY-IS with Mixture-based IS
# Runs everything and saves ONE list: mixture_results/mixture_results.rds
# Plots are made separately in mixture_output.R
###########################################

source("mixture_functions.R")

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

# list that will hold everything
results <- list()

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

laplace_gaussian_summary <- data.frame(
  sigma2 = laplace_gaussian_sigma2,
  mixture_ess = laplace_gaussian_output$diagnostics$ess,
  mixture_relative_ess = laplace_gaussian_output$diagnostics$relative_ess,
  my_weight_ess = laplace_gaussian_my$diagnostics$ess,
  my_weight_relative_ess = laplace_gaussian_my$diagnostics$relative_ess
)

# store everything for this example in the results list
results$laplace_gaussian <- list(
  lambda = laplace_gaussian_lambda,
  delta_px = laplace_gaussian_delta_px,
  delta_my = laplace_gaussian_delta_my,
  pilot = laplace_gaussian_pilot,
  means = laplace_gaussian_means,
  ess = laplace_gaussian_ess,
  mean_ess = laplace_gaussian_mean_ess,
  sigma2 = laplace_gaussian_sigma2,
  output = laplace_gaussian_output,
  my = laplace_gaussian_my,
  mix_logs = laplace_gaussian_mix_logs,
  my_logs = laplace_gaussian_my_logs,
  summary = laplace_gaussian_summary
)

# 2. Laplace / Student t

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

laplace_student_summary <- data.frame(
  sigma2 = laplace_student_sigma2,
  mixture_ess = laplace_student_output$diagnostics$ess,
  mixture_relative_ess = laplace_student_output$diagnostics$relative_ess,
  my_weight_ess = laplace_student_my$diagnostics$ess,
  my_weight_relative_ess = laplace_student_my$diagnostics$relative_ess
)

# store everything for this example in the results list
results$laplace_student <- list(
  lambda = laplace_student_lambda,
  delta_px = laplace_student_delta_px,
  delta_my = laplace_student_delta_my,
  pilot = laplace_student_pilot,
  means = laplace_student_means,
  ess = laplace_student_ess,
  mean_ess = laplace_student_mean_ess,
  sigma2 = laplace_student_sigma2,
  output = laplace_student_output,
  my = laplace_student_my,
  mix_logs = laplace_student_mix_logs,
  my_logs = laplace_student_my_logs,
  summary = laplace_student_summary
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

exp4_gaussian_summary <- data.frame(
  sigma2 = exp4_gaussian_sigma2,
  mixture_ess = exp4_gaussian_output$diagnostics$ess,
  mixture_relative_ess = exp4_gaussian_output$diagnostics$relative_ess,
  my_weight_ess = exp4_gaussian_my$diagnostics$ess,
  my_weight_relative_ess = exp4_gaussian_my$diagnostics$relative_ess
)

# store everything for this example in the results list
results$exp4_gaussian <- list(
  lambda = exp4_gaussian_lambda,
  delta_px = exp4_gaussian_delta_px,
  delta_my = exp4_gaussian_delta_my,
  pilot = exp4_gaussian_pilot,
  means = exp4_gaussian_means,
  ess = exp4_gaussian_ess,
  mean_ess = exp4_gaussian_mean_ess,
  sigma2 = exp4_gaussian_sigma2,
  output = exp4_gaussian_output,
  my = exp4_gaussian_my,
  mix_logs = exp4_gaussian_mix_logs,
  my_logs = exp4_gaussian_my_logs,
  summary = exp4_gaussian_summary
)

# 4. Quartic / Student t

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

exp4_student_summary <- data.frame(
  sigma2 = exp4_student_sigma2,
  mixture_ess = exp4_student_output$diagnostics$ess,
  mixture_relative_ess = exp4_student_output$diagnostics$relative_ess,
  my_weight_ess = exp4_student_my$diagnostics$ess,
  my_weight_relative_ess = exp4_student_my$diagnostics$relative_ess
)

# store everything for this example in the results list
results$exp4_student <- list(
  lambda = exp4_student_lambda,
  delta_px = exp4_student_delta_px,
  delta_my = exp4_student_delta_my,
  pilot = exp4_student_pilot,
  means = exp4_student_means,
  ess = exp4_student_ess,
  mean_ess = exp4_student_mean_ess,
  sigma2 = exp4_student_sigma2,
  output = exp4_student_output,
  my = exp4_student_my,
  mix_logs = exp4_student_mix_logs,
  my_logs = exp4_student_my_logs,
  summary = exp4_student_summary
)

# Save everything as a single list

results$settings <- list(
  d = d, beta = beta, start = start, df = df,
  pilot_n = pilot_n, pilot_thin = pilot_thin,
  sigma2_vec = sigma2_vec, iter_sigma2 = iter_sigma2,
  reps = reps, samp = samp
)

saveRDS(results, file.path(out_dir, "mixture_results.rds"))
