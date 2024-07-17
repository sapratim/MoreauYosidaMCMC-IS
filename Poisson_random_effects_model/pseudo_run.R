source("Poisson_functions.R")
#library(SimTools)

eta_start <- log(rowMeans(data)+1)
mu_start <- mean(eta_start)
lambda <- 0.001
rep <- 1e5

ismala <- mymala(eta_start, mu_start, lambda, sigma_eta, iter = rep, delta = 0.0032, data)
pxmala <- px.mala(eta_start, mu_start, lambda, sigma_eta, iter = rep, delta = 0.00065, data)
isbark <- mybarker(eta_start, mu_start, lambda, sigma_eta, iter = rep, delta = 0.003, data)
pxbark <- px.barker(eta_start, mu_start, lambda, sigma_eta, iter = rep, delta = 0.0006, data)

output_chain <- list(ismala, pxmala, isbark, pxbark)
save(output_chain, file = "output_pseudo.Rdata")