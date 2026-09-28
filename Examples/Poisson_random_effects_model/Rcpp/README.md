# Poisson random effects: Rcpp implementation

All seven samplers run in C++: `mymala`, `px.mala`, `mybarker`, `px.barker`,
`barker`, `myhmc`, and `pxhmc`. `Poisson_functions.R` compiles the engine with
`Rcpp::sourceCpp` and exposes the original R function names and positional
arguments. The original files in the parent directory are unchanged.

## Two separate workflows

Run scripts from the **Poisson_random_effects_model project root**, as in the
NuclearNorm example. Required R packages are `Rcpp`, `mcmcse`, `foreach`, and
`doParallel`, with an R-compatible C++ compiler. Runtime comparisons also use
`rbenchmark`.

Original R workflow:

```r
source("Poisson_run.R")
source("single_run_poisson.R")
source("Poisson_output.R")
```

Rcpp workflow:

```r
source("Rcpp/Poisson_run.R")
source("Rcpp/single_run_poisson.R")
source("Rcpp/Poisson_output.R")
```

The Rcpp scripts save all `.Rdata` files under `Rcpp` and the five plots under
`Rcpp/plots`. Both workflows keep the same saved object names, the same
17-element replication result, and the same ordering of single-chain outputs.
The output script derives the replication count from the saved results and
creates its plot directory. It does not clear the user's workspace.

Settings remain at the top of the run scripts: one million stored states,
100 replications, and 50 workers for the replication experiment. The data,
seed `8024248`, starting state, sampler parameters, and all enabled samplers
are preserved. The single-chain script retains its original PxHMC step size
of 0.0022; the replication script retains 0.0025.

A million-by-51 sample matrix alone occupies 408 MB. HMC also allocates a
momentum matrix of that size to preserve the original RNG order, and covariance
estimation requires additional matrices. Adjust `num_cores` to available memory;
50 workers can require tens of GB. Chains are released between sampler families
in the new replication runner.

Each replication run creates a fresh compilation cache. PSOCK workers source
the compiled functions sequentially from that cache, avoiding serialization of
native pointers and concurrent cache writes. Sampling then runs in parallel.
Workers receive separate random streams using `clusterSetRNGStream`; results
are repeatable with the same worker configuration. Clusters close on success
or error. Progress is forwarded to the console and can interleave.

## Direct calls and compatibility

```r
source("Rcpp/Poisson_functions.R")
eta_start <- log(rowMeans(data) + 1)
result <- mymala(eta_start, mean(eta_start), lambda = .001, sigma = sigma_eta,
                 iter = 1000, delta = .0032, data = data, verbose = FALSE)
samples <- result[[1]]
weights <- importance_weights(result[[2]])
```

| Sampler | Return value |
| --- | --- |
| `mymala` | `list(samples, log_weights)` |
| `px.mala` | sample matrix |
| `mybarker`, `myhmc` | `list(samples, log_weights, acceptance)` |
| `px.barker`, `barker`, `pxhmc` | `list(samples, acceptance)` |

Row 1 is the starting point. `verbose = TRUE` is the default, with ten progress
updates and a final acceptance report. `iter = 1` returns just the initial state.
Returned acceptance values retain `accepted / iter` for compatibility with
NuclearNorm and the original functions; printed rates use the actual number of
proposals, `iter - 1`. With no proposals, the printed rate is N/A and the returned
value is zero. `prior_sd` is an optional final argument defaulting to the original
`c = 10`. Model dimensions come from the supplied data, and all target/gradient
calculations use the supplied `sigma` and data.

The helper signatures are retained. `proxfunc` accepts the old initial-guess
arguments but uses a deterministic initializer based on its center instead.
The low-level exports are `poisson_sample_cpp` and `poisson_prox_cpp`; ordinary
calls should use the validated R wrappers. There are no R callbacks during the
sampling loops. Covariance estimation remains in R using `mcmcse`.

## Correctness review

The original MY target and log-weight sign are correct:
`log_weight = log p(x) - log p_lambda(x)`. The Px methods correctly use the
original target for acceptance and the smoothed gradient for proposals.
The original self-normalized importance-sampling covariance calculation is also
correct; its delta-method formula is retained.

The new implementation addresses these issues:

- The original log target uses global `sigma_eta`, while proposal gradients use
  the argument `sigma`. Gradients and proximal solves also read global `data`,
  `I`, and `ni_s`, even when a sampler receives other data. These discrepancies
  disappear for the default experiment, but can change the intended algorithm
  for other inputs. The engine consistently uses its explicit inputs.
- `log1p(exp(z))` in the original Barker correction can overflow. Stable
  softplus and logistic formulas remove this avoidable overflow. The omitted
  `dimension * log(2)` density constant cancels in MH and is not a sampler bug.
- The original proximal solver has no iteration limit or damping. The new
  solver uses the full arrowhead Hessian, an O(I) Schur-complement solve,
  backtracking, and a 100-iteration cap with explicit failure. It checks the
  stationarity residual against `tol_nr` plus a floating-point rounding allowance.
  Its initializer depends only on the current point, keeping numerical HMC
  forces independent of chain history.
- The original `2:iter` loop is invalid for `iter = 1`. The new engine handles
  that case and validates counts, dimensions, positive scales and iteration counts.
- Original printed acceptance rates divide by the number of stored states,
  rather than proposals. Printing is corrected; returned values deliberately
  retain the existing convention described above.
- Directly exponentiating very negative log weights can produce all zeros.
  `importance_weights` subtracts their maximum first; common rescaling leaves
  both self-normalized estimates and weight efficiency unchanged.

The new implementation has the same intended kernels and R random-number draw
order, including the upfront column-major HMC momentum draws and the 5% chance
of one leapfrog step. Different floating-point arithmetic and proximal solves
can eventually make long seeded trajectories diverge.

## Work avoided, especially for weights

Data row sums and fixed model coefficients are computed once per sampler call.
The current gradient and target are cached across iterations, including
rejections. MALA and Barker perform only one new proximal solve per proposal.
HMC performs one per leapfrog position and reuses the final solve for acceptance
and weights. The Newton solver avoids constructing or factoring dense Hessians.

MY weights need no additional proximal solve. They are evaluated only at the
initial state and after acceptance, using cached proximal exponentials. Direct
target differences using `expm1` reduce cancellation; a rejection copies the
previous log weight. The remaining cost is one O(I) pass for each accepted MY
state. The original samplers already avoided recomputing weights on rejection;
that behavior is preserved. Symmetric Gaussian terms in Barker and normalizing
constants in MALA cancel algebraically from the MH correction.

## Validation and timing

```sh
Rscript Rcpp/test_poisson.R
Rscript Rcpp/test_poisson.R --parallel
Rscript Rcpp/test_poisson.R --workflow
Rscript Rcpp/comparison.R
```

Regression checks cover all seven seeded R/Rcpp comparisons, per-state weight
identities and rejected-state reuse, `iter = 1`, independent BFGS and stationarity
checks for the prox, extreme centers, nondefault data and sigma, stable Barker
corrections, weight scaling, covariance agreement, and invalid iteration counts.
The parallel check tests all seven samplers on two workers and repeatability
with independent streams. The workflow check runs the actual replication,
single-chain and plotting scripts with reduced settings in a temporary project,
checking all saved layouts and the five generated plot files.

Initial validation passed: seeded differences were at most about 1.1e-9 over
80 stored states; the end-to-end check used two 4,000-state replications and
1,000-state single chains. `mcmcse` issued covariance definiteness warnings for
these deliberately short, highly correlated chains and used its batch-means
fallback. This is a workflow check, not a convergence assessment. The full
million-state, 100-replication experiment has not been run.

`comparison.R` uses `rbenchmark::benchmark` with 2,000 stored states and three
repetitions per sampler by default. Compilation, progress printing and variance
estimation are excluded; each timed call resets the same seed. It writes
`Rcpp/runtime_comparison.csv`. `elapsed` is total seconds across repetitions;
`relative` compares implementations within each sampler pair. Timings depend on
the machine and run length and do not predict full experiment runtime including
covariance estimation.

Measured on this machine on 2026-09-28 (three repetitions of 2,000 states):

| Sampler | R total seconds | Rcpp total seconds | Speedup |
| --- | ---: | ---: | ---: |
| MYMALA | 5.415 | 0.044 | 123.1x |
| PxMALA | 3.648 | 0.043 | 84.8x |
| MYBarker | 7.070 | 0.066 | 107.1x |
| PxBarker | 5.252 | 0.065 | 80.8x |
| Barker | 1.596 | 0.042 | 38.0x |
| MYHMC | 19.604 | 0.272 | 72.1x |
| PxHMC | 25.022 | 0.255 | 98.1x |

These are short sampler-only timings, including MY weight calculation.
