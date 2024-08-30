
source("TF_functions.R")
iter_warmup <- 1e5
lamb_coeff <- 0.001
D_mat <- getD(k=1, n=1e2, x)   #  D matrix
delta_warmup <- .0008
warmup_chain <- px.mala(y, alpha_hat, sigma2_hat, k=1, grid=x, iter = iter_warmup, 
                        delta = delta_warmup, start = y)
warmup_end_iter <- warmup_chain[iter_warmup,]
save(warmup_end_iter, file = "warmup.Rdata")
