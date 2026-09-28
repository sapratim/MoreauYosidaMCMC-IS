# Run from the project root: Rscript Rcpp/test_poisson.R [--parallel]
original <- new.env()
invisible(capture.output(sys.source("Poisson_functions.R", original)))
compiled <- new.env()
sys.source("Rcpp/Poisson_functions.R", compiled)
stopifnot(identical(original$data, compiled$data))
x <- log(rowMeans(compiled$data) + 1)
mu0 <- mean(x)
args <- list(eta_start = x, mu_start = mu0, lambda = .001, sigma = 3,
             iter = 80L, delta = .0032, data = compiled$data)
samplers <- c("mymala", "px.mala", "mybarker", "px.barker", "barker", "myhmc", "pxhmc")
arguments <- function(name, iterations = 80L) {
  a <- args; a$iter <- iterations
  if (name == "barker") a$lambda <- NULL
  if (name %in% c("myhmc", "pxhmc")) {
    a$delta <- NULL; a$eps_hmc <- if (name == "myhmc") .06 else .0025; a$L <- 10L
  }
  a
}
chain <- function(result) if (is.matrix(result)) result else result[[1]]
for (name in samplers) {
  a <- arguments(name)
  set.seed(746)
  invisible(capture.output(ref <- do.call(original[[name]], a)))
  set.seed(746)
  got <- do.call(compiled[[name]], c(a, list(verbose = FALSE)))
  error <- max(abs(unlist(ref) - unlist(got)))
  stopifnot(error < 2e-6)
  s <- chain(got)
  stopifnot(identical(dim(s), c(80L, 51L)), all(s[1, ] == c(x, mu0)))
  if (startsWith(name, "my")) {
    expected <- apply(s, 1, function(z) {
      p <- compiled$proxfunc(z[1:50], z[51], .001, x, mu0, 3)
      compiled$log_p(z[1:50], z[51], compiled$data) -
        compiled$log_plam(p[1:50], p[51], compiled$data, .001, z)
    })
    stopifnot(max(abs(expected - got[[2]])) < 1e-7, max(got[[2]]) <= 1e-8)
    repeated <- which(rowSums(abs(s[-1, ] - s[-nrow(s), ])) == 0) + 1L
    stopifnot(length(repeated) > 0, all(got[[2]][repeated] == got[[2]][repeated - 1L]))
  }
  one <- do.call(compiled[[name]], c(arguments(name, 1L), list(verbose = FALSE)))
  stopifnot(nrow(chain(one)) == 1L)
  if (is.list(one) && name != "mymala") stopifnot(tail(one, 1)[[1]] == 0)
  cat(name, "seeded maximum difference:", format(error), "\n")
}

# Independent proximal objective, finite-difference gradient, and BFGS oracle.
tiny <- matrix(c(0, 2, 10, 1, 4, 8), 3, 2)
z <- c(-2, .8, 3, .3); sigma <- 1.7; prior <- 2.2; lambda <- .2
objective <- function(p) sum((p[1:3] - p[4])^2)/(2*sigma^2) + p[4]^2/(2*prior^2) +
  2*sum(exp(p[1:3])) - sum(rowSums(tiny)*p[1:3]) + sum((p-z)^2)/(2*lambda)
p <- compiled$poisson_prox_cpp(z, tiny, lambda, sigma, prior)
opt <- optim(z, objective, method = "BFGS", control = list(reltol = 1e-13, maxit = 1000))
stopifnot(max(abs(p - opt$par)) < 2e-5)
g <- compiled$grad_logp(p[4], sigma, p[1:3], tiny, prior) - (p-z)/lambda
stopifnot(max(abs(g)) < 1e-7)
for (j in 1:4) {
  h <- rep(0,4); h[j] <- 1e-5
  stopifnot(abs((objective(p+h)-objective(p-h))/2e-5) < 1e-6)
}
# Extreme but finite centers and heterogeneous count scales.
for (center in list(c(-100, 100, 20, -50), c(1000, -1000, 40, 200))) {
  p <- compiled$poisson_prox_cpp(center, tiny, .001, sigma, prior)
  residual <- compiled$grad_logp(p[4], sigma, p[1:3], tiny, prior) - (p-center)/.001
  stopifnot(all(is.finite(p)), max(abs(residual)) < 1e-6)
}
# Explicit data and sigma must govern every sampler, regardless of globals.
for (name in samplers) {
  a <- arguments(name, 15L); a$eta_start <- z[1:3]; a$mu_start <- z[4]
  a$data <- tiny; a$sigma <- sigma; a$prior_sd <- prior; a$verbose <- FALSE
  if (!is.null(a$lambda)) a$lambda <- lambda
  set.seed(65)
  got <- do.call(compiled[[name]], a)
  stopifnot(identical(dim(chain(got)), c(15L, 4L)), all(is.finite(chain(got))))
  # Changing unrelated global model settings must not alter explicit-input calls.
  compiled$data <- matrix(999, 1, 1); compiled$I <- 1; compiled$ni_s <- 1
  compiled$sigma_eta <- 99
  set.seed(65)
  again <- do.call(compiled[[name]], a)
  stopifnot(identical(got, again))
  compiled$data <- original$data; compiled$I <- original$I
  compiled$ni_s <- original$ni_s; compiled$sigma_eta <- original$sigma_eta
  if (startsWith(name, "my")) {
    expected <- apply(chain(got), 1, function(state) {
      p <- compiled$poisson_prox_cpp(state, tiny, lambda, sigma, prior)
      compiled$log_p(state[1:3], state[4], tiny, sigma, prior) -
        compiled$log_plam(p[1:3], p[4], tiny, lambda, state, sigma, prior)
    })
    stopifnot(max(abs(expected - got[[2]])) < 1e-8)
  }
}
stopifnot(is.finite(compiled$log_bark.dens(0, 1, -1000, 1)))
stopifnot(all.equal(compiled$importance_weights(c(-10000, -10001, -Inf)), c(1, exp(-1), 0)))
# The original delta-method covariance is correct; rescaling preserves it.
set.seed(22); s <- matrix(rnorm(4000), 1000, 4); w <- exp(rnorm(1000))
stopifnot(isTRUE(all.equal(original$asymp_covmat_fn(s, w), compiled$asymp_covmat_fn(s, w), tolerance=1e-10)))
stopifnot(isTRUE(all.equal(compiled$asymp_covmat_fn(s, w), compiled$asymp_covmat_fn(s, 1e100*w), tolerance=1e-10)))
for (bad in list(0, -2, 1.5, NA_real_)) {
  a <- arguments("mymala"); a$iter <- bad; a$verbose <- FALSE
  stopifnot(inherits(try(do.call(compiled$mymala, a), silent=TRUE), "try-error"))
}
for (file in c("Poisson_run.R", "Poisson_output.R", "single_run_poisson.R")) parse(file.path("Rcpp",file))

if ("--parallel" %in% commandArgs(TRUE)) {
  cache <- tempfile("poisson_test_cache_")
  options(poisson.rcpp.cache_dir = cache)
  sys.source("Rcpp/Poisson_functions.R", compiled)
  cl <- parallel::makePSOCKcluster(2)
  tryCatch({
    for (i in seq_along(cl)) parallel::clusterCall(cl[i], function(root, cache) {
      setwd(root); options(poisson.rcpp.cache_dir = cache)
      source("Rcpp/Poisson_functions.R"); NULL
    }, getwd(), cache)
    doParallel::registerDoParallel(cl)
    run <- function() foreach::foreach(k = 1:2, .noexport = samplers) %dopar% {
      x <- log(rowMeans(data)+1)
      lapply(c("mymala", "px.mala", "mybarker", "px.barker", "barker", "myhmc", "pxhmc"), function(n) {
        a <- list(eta_start=x, mu_start=mean(x), sigma=3, iter=10L, data=data, verbose=FALSE)
        if (n != "barker") a$lambda <- .001
        if (n %in% c("myhmc", "pxhmc")) { a$eps_hmc <- .002; a$L <- 2L } else a$delta <- .001
        do.call(get(n), a)
      })
    }
    parallel::clusterSetRNGStream(cl, 42); first <- run()
    parallel::clusterSetRNGStream(cl, 42); second <- run()
    stopifnot(identical(first, second), !identical(first[[1]], first[[2]]))
  }, finally = { parallel::stopCluster(cl); foreach::registerDoSEQ() })
  cat("Parallel reproducibility passed.\n")
}
if ("--workflow" %in% commandArgs(TRUE)) {
  # Exercise the actual scripts with small settings in a disposable project.
  # Production settings and user outputs are never overwritten.
  root <- getwd(); sandbox <- tempfile("poisson_workflow_")
  dir.create(file.path(sandbox, "Rcpp"), recursive = TRUE)
  files <- c("Poisson_functions.R", "Poisson_functions.cpp", "Poisson_run.R",
             "Poisson_output.R", "single_run_poisson.R")
  file.copy(file.path(root, "Rcpp", files), file.path(sandbox, "Rcpp"))
  run_small <- function(file, settings) {
    env <- new.env(parent = globalenv())
    for (expr in parse(file)) {
      if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
          is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% names(settings)) {
        expr[[3]] <- settings[[as.character(expr[[2]])]]
      }
      eval(expr, env)
    }
    env
  }
  tryCatch({
    setwd(sandbox)
    run_small("Rcpp/Poisson_run.R", list(iter = 4000L, reps = 2L, num_cores = 2L))
    run_small("Rcpp/single_run_poisson.R", list(rep = 1000L))
    run_small("Rcpp/Poisson_output.R", list())
    saved <- new.env(); load("Rcpp/output_poisson.Rdata", saved)
    stopifnot(length(saved$output_poisson) == 2L,
              all(lengths(saved$output_poisson) == 17L))
    pdfs <- list.files("Rcpp/plots", pattern = "[.]pdf$", full.names = TRUE)
    stopifnot(length(pdfs) == 5L, all(file.info(pdfs)$size > 1000),
              !file.exists("output_poisson.Rdata"))
  }, finally = setwd(root))
  cat("Replication, single-chain and plotting workflows passed.\n")
}
cat("All Poisson checks passed.\n")
