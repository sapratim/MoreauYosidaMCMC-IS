
#####################################################################
################# Exponential output visualisation ##################
#####################################################################
set.seed(100)
library(mcmcse)
source("exp_d4_functions.R")
# These files must contain replicated results and the saved lamb_d* grids.
load("output_d1.Rdata")
load("output_d5.Rdata")
load("output_d10.Rdata")
load("output_d20.Rdata")
load("output_d50.Rdata")
load("output_d100.Rdata")

# Within EACH replication, choose lambda giving relative IS ESS closest to 0.50.
target_ess <- 0.40
dimensions <- c(1, 5, 10, 20, 50, 100)
mean_vec <- se_vec <- ci_vec <- numeric(length(dimensions))
selection_details <- vector("list", length(dimensions))

for (i in seq_along(dimensions)) {
  d <- dimensions[i]
  output <- get(paste0("output_exp_d", d))
  lambda_grid <- get(paste0("lamb_d", d))
  if (length(output) < 2L || !is.list(output[[1]][[1]])) {
    stop("Run the revised exp_run_d4.R first: replicated results are required.")
  }
  lambda_index <- vapply(output, function(rep_output) {
    relative_ess <- vapply(rep_output, function(l) l[[4]], numeric(1))
    which.min(abs(relative_ess - target_ess))
  }, integer(1))
  
  # Average efficiency over all components at each replication's selected lambda.
  eff_reps <- vapply(seq_along(output), function(r)
    mean(output[[r]][[lambda_index[r]]][[3]]), numeric(1))
  ess_reps <- vapply(seq_along(output), function(r)
    output[[r]][[lambda_index[r]]][[4]], numeric(1))
  selection_details[[i]] <- data.frame(
    dimension = d, replication = seq_along(output),
    selected_lambda = lambda_grid[lambda_index],
    relative_IS_ESS = ess_reps, relative_efficiency = eff_reps)
  
  mean_vec[i] <- mean(eff_reps)
  se_vec[i] <- sd(eff_reps) / sqrt(length(eff_reps))
  ci_vec[i] <- qt(0.975, df = length(eff_reps) - 1) * se_vec[i]
}

# Inspect the lambda and achieved n_e/n for every replication.
selection_details <- do.call(rbind, selection_details)
print(selection_details)
print(data.frame(dimension = dimensions,
                 mean_efficiency = mean_vec,
                 lower_95 = mean_vec - ci_vec,
                 upper_95 = mean_vec + ci_vec))
# These intervals describe the replication-specific lambda-selection procedure.

######## d vs Relative efficiency ########

dim <- dimensions
dir.create("plots", showWarnings = FALSE)

pdf(file = "plots/exp_dim_vs_re.pdf", width = 8, height = 6)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(dim, mean_vec, type = 'o', main = "Super Gaussian",
     ylim = range(mean_vec - ci_vec, mean_vec + ci_vec), xlab = "dimension", ylab = "Relative efficiency",
     cex.lab = 1.8,  cex.axis = 1.5,  cex.main = 1.8)
segments(dim, mean_vec - ci_vec,
         dim, mean_vec + ci_vec, lwd = 3, col = "blue")
dev.off()