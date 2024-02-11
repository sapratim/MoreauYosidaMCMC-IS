
source("nuclear_norm_functions.R")
load("start_value.Rdata")
iter <- 1e4 
lamb_coeff <- 0.00001
sigma2_hat <- 0.01
alpha_hat <- 1.15/sigma2_hat
step_ismala <- 0.00008
step_pxmala <- 0.00008

result_is <- mymala(y=y, alpha = alpha_hat, sigma2 = sigma2_hat,
                      iter = iter, delta = step_ismala)
result_pxm <- px.mala(y=y, alpha = alpha_hat, sigma2 = sigma2_hat,
                      iter = iter, delta = step_pxmala)
result_ishmc <- myhmc(y=y, alpha = alpha_hat, sigma2 = sigma2_hat, 
                      iter = iter, eps_hmc = 0.0047, L=10)
result_pxhmc <- pxhmc(y=y, alpha = alpha_hat, sigma2 = sigma2_hat, 
                      iter = iter, eps_hmc = 0.004, L=10)

## For ISMALA
ismala_chain <- result_is[[1]]
ismala_weights <- result_is[[2]]
hist(ismala_weights)
traceplot(ismala_chain[, c(11:20)])
traceplot(ismala_weights)

## For PxMALA
pxmala_chain <- result_pxm[[1]]
pxmala_weights <- result_pxm[[2]]
hist(pxmala_weights)
traceplot(pxmala_chain[, c(11:20)])
traceplot(pxmala_weights)

## For ISHMC
ishmc_chain <- result_ishmc[[1]]
ishmc_weights <- result_ishmc[[2]]
traceplot(ishmc_chain[, c(11:20)])
traceplot(ishmc_weights)
hist(ishmc_weights)

## For ISHMC
pxhmc_chain <- result_ishmc[[1]]
pxhmc_weights <- result_ishmc[[2]]
traceplot(pxhmc_chain[, c(91:100)])
traceplot(pxhmc_weights)
hist(pxhmc_weights)

############# Posterior mean ISHMC

weight_mat <- matrix(0, nrow = 1e4, ncol = length(y))
for (i in 1:1e4) {
  weight_mat[i,] <- ishmc_chain[i,]*exp(ishmc_weights[i])
}
num_sum <- apply(weight_mat, 2, sum)
weights_sum <- sum(exp(weights))
post_mean <- num_sum/weights_sum    #### Mean value

############# Posterior mode ISHMC
true_log_target <- function(image_vec)
{
  sigma2 <- 0.01
  alpha <- 1.15/sigma2
  n_norm <- nucl_norm(image_vec)
  dens_val <- alpha*n_norm + sum((y-image_vec)^2)/(2*sigma2)
  return(-dens_val)
}
post_val <- numeric(length = iter)
for (i in 1:iter) {
  post_val[i] <- true_log_target(ishmc_chain[i,])
}
mode_index <- which.max(post_val)

post_mode1 <- ishmc_chain[mode_index,]   ### Mode value 

############ Image visualisation
pdf(file = "sample_plot.pdf", width = 6, height = 4)
par(mfrow = c(1,2))
image(matrix(post_mode1, nrow = n, ncol = n), col = gray.colors(4, start = 0, end = 1), axes = FALSE)
image(checker, col = gray.colors(4, start = 0, end = 1), axes = FALSE)
dev.off()
