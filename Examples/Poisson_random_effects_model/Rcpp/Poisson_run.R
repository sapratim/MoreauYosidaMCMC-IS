######## Poisson random effects run ##########

# Run from the project root. Use a fresh shared compilation cache per run.
cache_dir <- tempfile("poisson_rcpp_")
options(poisson.rcpp.cache_dir = cache_dir)
source("Rcpp/Poisson_functions.R")

lambda <- 0.001
eta_start <- log(rowMeans(data)+1)
mu_start <- mean(eta_start)
iter <- 1e6
delta_mym <- 0.0032
delta_pxm <- 0.00065
delta_mybark <- 0.003
delta_pxbark <- 0.0006
delta_bark <- 0.0012

output_poisson <- list()
parallel::detectCores()
num_cores <- 50
reps <- 100

# Load compiled functions in each worker; never serialize Rcpp pointers.
# Sequential cache loads avoid simultaneous sourceCpp cache writes.
cluster <- parallel::makePSOCKcluster(num_cores, outfile = "")
output_poisson <- tryCatch({
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

output_poisson <- foreach(b = 1:reps, .packages = "mcmcse",
  .noexport = c("mymala", "px.mala", "mybarker", "px.barker", "barker", "myhmc", "pxhmc",
                "poisson_sample_cpp", "asymp_covmat_fn", "importance_weights")) %dopar% {
  
################  MALA  ##################
  
  time_ism <- system.time(ismala <- mymala(eta_start, mu_start, lambda, sigma_eta, 
                                          iter = iter, delta = delta_mym, data))
  time_pxm <- system.time(pxmala <- px.mala(eta_start, mu_start, lambda, sigma_eta,
                                           iter = iter, delta = delta_pxm, data))
  
  time_ism <- time_ism["elapsed"]
  time_pxm <- time_pxm["elapsed"]
  
  mala_chain <- ismala[[1]]
  weights_ism <- importance_weights(ismala[[2]])
  n_eff_mala <- (mean(weights_ism)^2)/mean(weights_ism^2)
  
  # Asymptotic covariance matrix
  asymp_covmat_ism <- asymp_covmat_fn(mala_chain, weights_ism) 
  asymp_covmat_pxm <- mcse.multi(pxmala)$cov   
  
rm(ismala, pxmala, mala_chain, weights_ism)

################  Barker  ##################  
  
  time_isb <- system.time(isbark <- mybarker(eta_start, mu_start, lambda, sigma_eta, 
                                             iter = iter, delta = delta_mybark, data))
  time_pxb <- system.time(pxbark <- px.barker(eta_start, mu_start, lambda, sigma_eta,
                                              iter = iter, delta = delta_pxbark, data))
  time_trueb <- system.time(true_bark <- barker(eta_start, mu_start, sigma_eta, 
                                                iter = iter, delta = delta_bark, data))
  
  time_isb <- time_isb["elapsed"]
  time_pxb <- time_pxb["elapsed"]
  time_trueb <- time_trueb["elapsed"]
  
  bark_chain <- isbark[[1]]
  weights_isb <- importance_weights(isbark[[2]])
  n_eff_bark <- (mean(weights_isb)^2)/mean(weights_isb^2)
  
  # Asymptotic covariance matrix
  asymp_covmat_isb <- asymp_covmat_fn(bark_chain, weights_isb) 
  asymp_covmat_pxb <- mcse.multi(pxbark[[1]])$cov   
  asymp_covmat_trubark <- mcse.multi(true_bark[[1]])$cov
  
rm(isbark, pxbark, true_bark, bark_chain, weights_isb)

################  HMC  ##################    
  
  time_ish <- system.time(my.hmc <- myhmc(eta_start, mu_start,lambda, 
                                      sigma_eta, iter = iter, data, eps_hmc=0.06, L=10))
  time_pxh <- system.time(px.hmc <- pxhmc(eta_start, mu_start,lambda, 
                                      sigma_eta, iter = iter, data, eps_hmc=0.0025, L=10)) 
  
  time_ish <- time_ish["elapsed"]
  time_pxh <- time_pxh["elapsed"]
  
  hmc_chain <- my.hmc[[1]]
  weights_hmc <- importance_weights(my.hmc[[2]])
  n_eff_hmc <- (mean(weights_hmc)^2)/mean(weights_hmc^2)
  
  # Asymptotic covariance matrix
  asymp_covmat_ishmc <- asymp_covmat_fn(hmc_chain, weights_hmc) 
  asymp_covmat_pxhmc <- mcse.multi(px.hmc[[1]])$cov   # PxMALA asymptotic variance
  
  list(asymp_covmat_ism, asymp_covmat_pxm, asymp_covmat_isb, asymp_covmat_pxb, 
       asymp_covmat_trubark, asymp_covmat_ishmc, asymp_covmat_pxhmc, time_ism, 
       time_pxm, time_isb, time_pxb, time_trueb, time_ish, time_pxh, n_eff_mala,
                                                               n_eff_bark, n_eff_hmc)
}

output_poisson
}, finally = {
  parallel::stopCluster(cluster)
  foreach::registerDoSEQ()
})
save(output_poisson, file = "Rcpp/output_poisson.Rdata")
