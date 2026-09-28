# Run from the project root: source("Rcpp/comparison.R")
# Compilation and variance estimation are excluded from these timings.
original <- new.env()
invisible(capture.output(sys.source("Poisson_functions.R", original)))
compiled <- new.env()
sys.source("Rcpp/Poisson_functions.R", compiled)
# Silence original progress without changing its sampler calculations.
original$print <- function(x, ...) invisible(x)
original$cat <- function(...) invisible(list(...))
iter <- 2000
reps <- 3
eta_start <- log(rowMeans(compiled$data) + 1)
mu_start <- mean(eta_start)
settings <- list(mymala = .0032, px.mala = .00065, mybarker = .003,
                 px.barker = .0006, barker = .0012, myhmc = .06, pxhmc = .0025)
results <- lapply(names(settings), function(name) {
  args <- list(eta_start = eta_start, mu_start = mu_start, sigma = 3,
               iter = iter, data = compiled$data)
  if (name != "barker") args$lambda <- .001
  if (name %in% c("myhmc", "pxhmc")) {
    args$eps_hmc <- settings[[name]]; args$L <- 10L
  } else args$delta <- settings[[name]]
  r_args <- args
  cpp_args <- c(args, list(verbose = FALSE))
  run_r <- function() { set.seed(123); do.call(original[[name]], r_args) }
  run_cpp <- function() { set.seed(123); do.call(compiled[[name]], cpp_args) }
  invisible(run_cpp())
  b <- rbenchmark::benchmark(R = run_r(), Rcpp = run_cpp(), replications = reps,
                            columns = c("test", "replications", "elapsed", "relative"))
  data.frame(sampler = name, b, row.names = NULL)
})
runtime_comparison <- do.call(rbind, results)
print(runtime_comparison, row.names = FALSE)
write.csv(runtime_comparison, "Rcpp/runtime_comparison.csv", row.names = FALSE)
