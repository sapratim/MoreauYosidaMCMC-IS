
source("Poisson_functions.R")
eta_start <- log(rowMeans(data)+1)
mu_start <- mean(eta_start)
lambda <- 0.001
rep <- 1e5

ismala <- mymala(eta_start, mu_start, lambda, sigma_eta, iter = rep, delta = 0.0032, data)
pxmala <- px.mala(eta_start, mu_start, lambda, sigma_eta, iter = rep, delta = 0.00065, data)
isbark <- mybarker(eta_start, mu_start, lambda, sigma_eta, iter = rep, delta = 0.003, data)
pxbark <- px.barker(eta_start, mu_start, lambda, sigma_eta, iter = rep, delta = 0.0006, data)
true_bark <- barker(eta_start, mu_start, sigma_eta, iter = rep, delta = 0.0012, data)
my.hmc <- myhmc(eta_start, mu_start,lambda, sigma_eta, iter = rep, data, eps_hmc=0.06, L=10)
px.hmc <- pxhmc(eta_start, mu_start,lambda, sigma_eta, iter = rep, data, eps_hmc=0.002, L=10) 

output_chain_mala <- list(ismala[[1]], pxmala)
output_chain_bark <- list(isbark[[1]], pxbark[[1]], true_bark[[1]])
output_chain_hmc <- list(my.hmc[[1]], px.hmc[[1]])

save(output_chain_mala, file = "output_single_mala.Rdata")
save(output_chain_bark, file = "output_single_bark.Rdata")
save(output_chain_hmc, file = "output_single_hmc.Rdata")

