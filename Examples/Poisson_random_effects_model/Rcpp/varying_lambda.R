######## Poisson random effects run ##########

# Run from the project root. Use a fresh shared compilation cache per run.
cache_dir <- tempfile("poisson_rcpp_")
options(poisson.rcpp.cache_dir = cache_dir)
source("Rcpp/Poisson_functions.R")

lambda <- seq(0.0001, 0.01, length = 10)
eta_start <- log(rowMeans(data)+1)
mu_start <- mean(eta_start)
iter <- 1e4

delta_mym_seq <- c(0.00135, 0.0038, 0.0055, 0.0075, 0.009, 0.012, 0.013, 0.0135, 0.015, 0.017)
delta_pxm_seq <- c(0.0011, 0.00062, 0.0005, 0.00046, 0.00042, 0.00035, 0.00037, 0.00032, 0.0003, 0.00027)
delta_mybark_seq <- c(0.0015, 0.0032, 0.0045, 0.0058, 0.0068, 0.008, 0.0088, 0.01, 0.011, 0.0122)
delta_pxbark_seq <- c(0.0012, 0.00065, 0.0005, 0.00035, 0.00034, 0.00038, 0.00038, 0.0003, 0.0003, 0.00028)
delta_bark <- 0.0012

output_lambda <- list()
parallel::detectCores()
num_cores <- 20
reps <- 20

# Load compiled functions in each worker; never serialize Rcpp pointers.
# Sequential cache loads avoid simultaneous sourceCpp cache writes.
cluster <- parallel::makePSOCKcluster(num_cores, outfile = "")
output_lambda <- tryCatch({
  for (worker in seq_along(cluster)) {
    parallel::clusterCall(cluster[worker], function(project_dir, cache) {
      setwd(project_dir)
      options(poisson.rcpp.cache_dir = cache)
      source("Rcpp/Poisson_functions.R")
      NULL
    }, getwd(), cache_dir)
  }
  doParallel::registerDoParallel(cluster)
  parallel::clusterSetRNGStream(cluster, iseed = 20260928)
  
        output_lambda <- foreach(b = 1:reps, .packages = "mcmcse",
                            .noexport = c("mymala", "px.mala", "mybarker", "px.barker", "barker", "myhmc", "pxhmc",
                                          "poisson_sample_cpp", "asymp_covmat_fn", "importance_weights")) %do% {
                                            
        mala_lambda_var <- matrix(0, nrow = length(lambda), ncol = length(eta_start)+1)
        bark_lambda_var <- matrix(0, nrow = length(lambda), ncol = length(eta_start)+1)
        hmc_lambda_var <- matrix(0, nrow = length(lambda), ncol = length(eta_start)+1)
                                                                                
        time_trueb <- system.time(true_bark <- barker(eta_start, mu_start, sigma_eta, 
                                               iter = iter, delta = delta_bark, data))                                    
        time_trueb <- time_trueb["elapsed"]
        
        asymp_var_trubark <- diag(mcse.multi(true_bark[[1]])$cov)                                    
        rm(true_bark)                                    
                        
        for (i in 1:length(lambda)) {
          
             ################  MALA  ##################
        time_ism <- system.time(ismala <- mymala(eta_start, mu_start, lambda[i], sigma_eta, 
                                    iter = iter, delta = delta_mym_seq[i], data)) 
        time_pxm <- system.time(pxmala <- px.mala(eta_start, mu_start, lambda[i], sigma_eta,
                                 iter = iter, delta = delta_pxm_seq[i], data))
                                            
        time_ism <- time_ism["elapsed"]
        time_pxm <- time_pxm["elapsed"]
                                            
        mala_chain <- ismala[[1]]
        weights_ism <- importance_weights(ismala[[2]])
        n_eff_mala <- (mean(weights_ism)^2)/mean(weights_ism^2)
                                            
        # Asymptotic covariance matrix
       asymp_var_ism <- diag(asymp_covmat_fn(mala_chain, weights_ism)) 
       asymp_var_pxm <- diag(mcse.multi(pxmala)$cov)   
       
       mala_lambda_var[i,] <- (asymp_var_pxm*time_pxm)/(asymp_var_ism*time_ism)
       bark_lambda_var[i,] <- (asymp_var_trubark*time_trueb)/(asymp_var_ism*time_ism)
                                            
       rm(ismala, pxmala, mala_chain, weights_ism)
                                            
       ################  Barker  ##################  
                                            
     # time_isb <- system.time(isbark <- mybarker(eta_start, mu_start, lambda, sigma_eta, 
     #                                               iter = iter, delta = delta_mybark, data))
     # time_pxb <- system.time(pxbark <- px.barker(eta_start, mu_start, lambda, sigma_eta,
     #                                              iter = iter, delta = delta_pxbark, data))
     # 
     #                                        
     #  time_isb <- time_isb["elapsed"]
     #  time_pxb <- time_pxb["elapsed"]
     #  time_trueb <- time_trueb["elapsed"]
     #                                        
     #    bark_chain <- isbark[[1]]
     #    weights_isb <- importance_weights(isbark[[2]])
     #    n_eff_bark <- (mean(weights_isb)^2)/mean(weights_isb^2)
     #                                        
     #                                        # Asymptotic covariance matrix
     #  asymp_covmat_isb <- asymp_covmat_fn(bark_chain, weights_isb) 
     #  asymp_covmat_pxb <- mcse.multi(pxbark[[1]])$cov   
     #  asymp_covmat_trubark <- mcse.multi(true_bark[[1]])$cov
      
      
         ################  HMC  ##################    
     eps_is <- c(0.034, 0.058, 0.075, 0.08, 0.102, 0.102, 0.108, 0.12, 0.125, 0.132)
     eps_px <- c(0.03, 0.002, 0.0017, 0.0016, 0.0015, 0.00148, 0.00142, 0.0014, 0.0014, 0.00138)  
     
     time_ish <- system.time(my.hmc <- myhmc(eta_start, mu_start,lambda[i], 
                             sigma_eta, iter = iter, data, eps_hmc= eps_is[i], L=10))
      time_pxh <- system.time(px.hmc <- pxhmc(eta_start, mu_start,lambda[i], 
                            sigma_eta, iter = iter, data, eps_hmc= eps_px[i], L=10)) 
                                            
      time_ish <- time_ish["elapsed"]
      time_pxh <- time_pxh["elapsed"]
                                            
       hmc_chain <- my.hmc[[1]]
       weights_hmc <- importance_weights(my.hmc[[2]])
       n_eff_hmc <- (mean(weights_hmc)^2)/mean(weights_hmc^2)
                                            
       # Asymptotic covariance matrix
       asymp_var_ish <- diag(asymp_covmat_fn(hmc_chain, weights_hmc)) 
       asymp_var_pxh <- diag(mcse.multi(px.hmc[[1]])$cov)   # PxMALA asymptotic variance
       
       hmc_lambda_var[i,] <- (asymp_var_pxh*time_pxh)/(asymp_var_ish*time_ish)
       
        }                                         
  list(mala_lambda_var, bark_lambda_var, hmc_lambda_var, time_ism, time_pxm, 
         time_trueb, time_ish, time_pxh, n_eff_mala, n_eff_hmc)
  }
  
  output_lambda
}, finally = {
  parallel::stopCluster(cluster)
  foreach::registerDoSEQ()
})
save(output_lambda, file = "Rcpp/output_variable_lambda.Rdata")
