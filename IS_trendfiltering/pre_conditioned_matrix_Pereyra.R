source("IS_trendf_Pereyra.R")

iter <- 1e5

lamb_coeff <- 0.0005
D_mat <- getD(k=1, n=1e2, x)   #  D matrix
delta <- .0005
prec_mat_chain <- mymala_cov_fn(y, alpha_hat, sigma2_hat, k=1, 
                                            grid = x, iter, delta = delta)
markov_chain <- prec_mat_chain[[1]]
covmat <- prec_mat_chain[[2]]
save(markov_chain, file = "MC_pcm.Rdata")
save(covmat, file = "covmat.Rdata")
