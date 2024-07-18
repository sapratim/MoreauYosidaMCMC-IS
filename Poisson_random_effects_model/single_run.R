
source("Poisson_functions.R")
eta_start <- log(rowMeans(data)+1)
mu_start <- mean(eta_start)
lambda <- 0.001
rep <- 1e5

ismala <- mymala(eta_start, mu_start, lambda, sigma_eta, iter = rep, delta = 0.0032, data)
pxmala <- px.mala(eta_start, mu_start, lambda, sigma_eta, iter = rep, delta = 0.00065, data)
isbark <- mybarker(eta_start, mu_start, lambda, sigma_eta, iter = rep, delta = 0.003, data)
pxbark <- px.barker(eta_start, mu_start, lambda, sigma_eta, iter = rep, delta = 0.0006, data)
my.hmc <- myhmc(eta_start, mu_start,lambda, sigma_eta, iter = rep, data, eps_hmc=0.06, L=10)
px.hmc <- pxhmc(eta_start, mu_start,lambda, sigma_eta, iter = rep, data, eps_hmc=0.002, L=10) 

output_chain <- list(ismala[[1]], pxmala, isbark[[1]], pxbark[[1]], my.hmc[[1]], px.hmc[[1]])
save(output_chain, file = "output_single.Rdata")