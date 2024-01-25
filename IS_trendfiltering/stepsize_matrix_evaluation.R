
source("IS_trendf_functions_Pereyra.R")
iter_pcm <- 1e5
lamb_coeff <- 0.001
D_mat <- getD(k=1, n=1e2, x)   #  D matrix
delta_pcm <- .0015

prec_mat_chain <- mymala_cov_fn(y, alpha_hat, sigma2_hat, k=1, 
                                grid = x, iter_pcm, delta = delta_pcm)
markov_chain <- prec_mat_chain
covmat <- cov(markov_chain[-c(1:5e4),])
save(markov_chain, file = "MC_pcm.Rdata")
save(covmat, file = "covmat.Rdata")
