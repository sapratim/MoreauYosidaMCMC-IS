######## Poisson random effects run ##########

source("Poisson_functions.R")

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
doParallel::registerDoParallel(cores = num_cores)
reps <- 100

output_poisson <- foreach(b = 1:reps) %dopar% {
  
################  MALA  ##################
  
  time_ism <- system.time(ismala <- mymala(eta_start, mu_start, lambda, sigma_eta, 
                                          iter = iter, delta = delta_mym, data))
  time_pxm <- system.time(pxmala <- px.mala(eta_start, mu_start, lambda, sigma_eta,
                                           iter = iter, delta = delta_pxm, data))
  
  time_ism <- time_ism["elapsed"]
  time_pxm <- time_pxm["elapsed"]
  
  mala_chain <- matrix(unlist(ismala[[1]]), nrow = iter, ncol = I+1)
  weights_ism <- exp(as.numeric(unlist(ismala[[2]])))
  n_eff_mala <- (mean(weights_ism)^2)/mean(weights_ism^2)
  
  # Asymptotic covariance matrix
  asymp_covmat_ism <- asymp_covmat_fn(mala_chain, weights_ism) 
  asymp_covmat_pxm <- mcse.multi(pxmala)$cov   
  
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
  
  bark_chain <- matrix(unlist(isbark[[1]]), nrow = iter, ncol = I+1)
  weights_isb <- exp(as.numeric(unlist(isbark[[2]])))
  n_eff_bark <- (mean(weights_isb)^2)/mean(weights_isb^2)
  
  # Asymptotic covariance matrix
  asymp_covmat_isb <- asymp_covmat_fn(bark_chain, weights_isb) 
  asymp_covmat_pxb <- mcse.multi(pxbark[[1]])$cov   
  asymp_covmat_trubark <- mcse.multi(true_bark[[1]])$cov
  
################  HMC  ##################    
  
  time_ish <- system.time(my.hmc <- myhmc(eta_start, mu_start,lambda, 
                                      sigma_eta, iter = iter, data, eps_hmc=0.06, L=10))
  time_pxh <- system.time(px.hmc <- pxhmc(eta_start, mu_start,lambda, 
                                      sigma_eta, iter = iter, data, eps_hmc=0.0025, L=10)) 
  
  time_ish <- time_ish["elapsed"]
  time_pxh <- time_pxh["elapsed"]
  
  hmc_chain <- matrix(unlist(my.hmc[[1]]), nrow = iter, ncol = I+1)
  weights_hmc <- exp(as.numeric(unlist(my.hmc[[2]])))
  n_eff_hmc <- (mean(weights_hmc)^2)/mean(weights_hmc^2)
  
  # Asymptotic covariance matrix
  asymp_covmat_ishmc <- asymp_covmat_fn(hmc_chain, weights_hmc) 
  asymp_covmat_pxhmc <- mcse.multi(px.hmc[[1]])$cov   # PxMALA asymptotic variance
  
  list(asymp_covmat_ism, asymp_covmat_pxm, asymp_covmat_isb, asymp_covmat_pxb, 
       asymp_covmat_trubark, asymp_covmat_ishmc, asymp_covmat_pxhmc, time_ism, 
       time_pxm, time_isb, time_pxb, time_trueb, time_ish, time_pxh, n_eff_mala,
                                                               n_eff_bark, n_eff_hmc)
}

save(output_poisson, file = "output_poisson.Rdata")
