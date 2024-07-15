
source("TF_functions.R")
load("pcm_last_iter.Rdata")

iter_mala <- 1e4
lamb_coeff <- 0.001
D_mat <- getD(k=1, n=1e2, x)   #  D matrix
delta_samp_is <- 0.0015
delta_samp_pxm <- 0.0008

mala.is <- mymala(y, alpha_hat, sigma2_hat, k=1, grid=x, iter = iter_mala, 
                  delta = delta_samp_is, start = pcm_last_iter)

pxmala.run <- px.mala(y, alpha_hat, sigma2_hat, k=1, grid=x, iter = iter_mala, 
                      delta = delta_samp_pxm, start = pcm_last_iter)

mala_chain <- mala.is[[1]]
weights_mala <- mala.is[[2]]
is_samp <- matrix(unlist(mala_chain), nrow = iter_mala, ncol = length(y))
is_wts <- as.numeric(unlist(weights))
wts_mean <- mean(exp(is_wts))
num <- is_samp*exp(is_wts)
sum_mat <- apply(num, 2, sum)
is_est <- sum_mat / sum(exp(is_wts))
input_mat <- cbind(num, exp(is_wts))  # input samples for mcse
Sigma_mat <- mcse.multi(input_mat)$cov  # estimated covariance matrix of the tuple
kappa_eta_mat <- cbind(diag(1/wts_mean, length(y)), -is_est/wts_mean) # derivative of kappa matrix

asymp_covmat_is <- (kappa_eta_mat %*% Sigma_mat) %*% t(kappa_eta_mat) # IS asymptotic variance

asymp_covmat_pxm <- mcse.multi(pxmala.run)$cov   # PxMALA asymptotic variance

rel_ess <- (det(asymp_covmat_pxm)/det(asymp_covmat_is))^(1/length(y))

load("pcm_last_iter.Rdata")

source("TF_functions.R")


iter_hmc <- 1e3
lamb_coeff <- 0.001
D_mat <- getD(k=1, n=1e2, x)   #  D matrix

  my.hmc <- myhmc(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter_hmc,
                  eps_hmc = 0.015, L = 100, start = pcm_last_iter)
  
  px.hmc <- pxhmc(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter_hmc,
                  eps_hmc = 0.0003, L = 100, start = pcm_last_iter) 
  
  hmc_chain <- my.hmc[[1]]
  weights <- my.hmc[[2]]
  pxhmc_chain <- px.hmc[[1]]
  is_samp <- matrix(unlist(hmc_chain), nrow = iter_hmc, ncol = length(y))
  is_wts <- as.numeric(unlist(weights))
  wts_mean <- mean(exp(is_wts))
  num <- is_samp*exp(is_wts)
  sum_mat <- apply(num, 2, sum)
  is_est <- sum_mat / sum(exp(is_wts))
  input_mat <- cbind(num, exp(is_wts))  # input samples for mcse
  Sigma_mat <- mcse.multi(input_mat)$cov  # estimated covariance matrix of the tuple
  kappa_eta_mat <- cbind(diag(1/wts_mean, length(y)), -is_est/wts_mean) # derivative of kappa matrix
  
  asymp_covmat_is <- (kappa_eta_mat %*% Sigma_mat) %*% t(kappa_eta_mat) # IS asymptotic variance
  
  asymp_covmat_pxhmc <- mcse.multi(pxhmc_chain)$cov   # PxMALA asymptotic variance
  
  rel_ess <- (det(asymp_covmat_pxhmc)/det(asymp_covmat_is))^(1/length(y))
  
  ##  Posterior mean
  
  weight_mat <- matrix(0, nrow = iter_hmc, ncol = length(y))
  for (i in 1:iter_hmc) {
    weight_mat[i,] <- hmc_chain[i,]*exp(weights[i])
  }
  num_sum <- apply(weight_mat, 2, sum)
  weights_sum <- sum(exp(weights))
  post_mean <- num_sum/weights_sum
  
  #  Quantile visualisation
  
  augm_mat <- cbind(hmc_chain,weights)
  
  upper_quant <- numeric(length = length(y))
  lower_quant <- numeric(length = length(y))
  post_med <- numeric(length = length(y))
  signif_level <- 0.025
  
  for (i in 1:length(y)) 
  {
    initial_mat <- quant(i, augm_mat)
    mat_sum <- apply(initial_mat, 2, sum)
    wts_prop <- initial_mat[,2]/mat_sum[2]
    final_mat <- cbind(initial_mat[1,], cumsum(wts_prop))
    lower_index <- min(which(final_mat[,2] >= signif_level))
    upper_index <- min(which(final_mat[,2] >= (1 - signif_level)))
    med_index <- min(which(final_mat[,2] >= 0.5))
    upper_quant[i] <- initial_mat[upper_index,1]
    lower_quant[i] <- initial_mat[lower_index,1]
    post_med[i] <- initial_mat[med_index,1]
  }
  
  acc_rate_is <- my.hmc[[3]]
  acc_rate_pxhmc <- px.hmc[[2]]
  list(post_mean, post_med, Sigma_mat, asymp_covmat_is, asymp_covmat_pxhmc, 
       upper_quant, lower_quant, rel_ess, acc_rate_is, acc_rate_pxhmc)

  
  #################   Barker
  
  iter_bark <- 1e4
  lamb_coeff <- 0.001
  D_mat <- getD(k=1, n=1e2, x)   #  D matrix
  delta_bark_my <- 0.06
  delta_bark_px <- 0.042
  
  
  pxbark <- px.barker(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter_bark,delta = delta_bark_px,
                     start = pcm_last_iter)
  mybark <- mybarker(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter_bark,delta = delta_bark_my,
                     start = pcm_last_iter)
  
  bark_chain <- mybark[[1]]
  weights <- mybark[[2]]

  
  
  
  
  
  
  
  
  
  
  #### MALA for pi-lambda without weights, for warm-up
  
  mymala_cov_fn <- function(y, alpha, sigma2, k, grid, iter, delta) # pre-conditioned mala
  {
    nvar <- length(y)
    samp.mym <- matrix(0, nrow = iter, ncol = nvar)
    lambda <- lamb_coeff
    
    # starting value computations
    beta_current <- y
    prox_val.curr <- prox_func(beta_current, lambda, alpha,sigma2, k, grid)
    targ_val.curr <- log_pilambda(prox_val.curr,beta_current,lambda,y,sigma2,alpha)
    samp.mym[1,] <- beta_current
    accept <- 0
    
    # mcmc
    for (i in 2:iter) 
    {
      
      # proposal
      prop.mean <- beta_current + (delta / 2)*grad_logpiLam(beta_current,lambda,y,sigma2,alpha,k,grid)
      beta_next <- rnorm(nvar, prop.mean, sd = sqrt(delta))
      
      # Evaluating proximal map and target (pi-lambda)
      prox_val.next <- prox_func(beta_next, lambda, alpha, sigma2, k, grid)
      targ_val.next <- log_pilambda(prox_val.next,beta_next,lambda,y,sigma2,alpha)
      
      other.mean <-  beta_next + (delta / 2)*grad_logpiLam(beta_next,lambda,y,sigma2,alpha,k,grid)
      q.next_to_curr <- sum(dnorm(beta_current, other.mean, sd = sqrt(delta), log = TRUE))
      q.curr_to_next <- sum(dnorm(beta_next, prop.mean, sd = sqrt(delta), log = TRUE))
      
      mh.ratio <- targ_val.next + q.next_to_curr - (targ_val.curr + q.curr_to_next)  # mh  ratio
      # print(mh.ratio)
      if(log(runif(1)) <= mh.ratio)
      {
        samp.mym[i,] <- beta_next
        prox_val.curr <- prox_val.next
        targ_val.curr <- targ_val.next
        accept <- accept + 1
      }
      else
      {
        samp.mym[i,] <- beta_current
      }
      beta_current <- samp.mym[i,]
      if(i %% 1000 == 0){
        j <- accept/iter
        print(cat(i, j))
      }
    }
    print(accept/iter)
    object <- samp.mym
    return(object)
  }

  
  
  bark.dens <- function(curr_point, prop_point, grad_curr_point)
  {
    rw_dens <- prod(dnorm(prop_point-curr_point, 0, 1))
    exp_term <- - (grad_curr_point*(prop_point-curr_point))
    denom <- prod(1 + exp(exp_term))
    dens_val <- rw_dens/denom
    return(dens_val)
  }
  
  
  
  mala_chain <- matrix(unlist(mala.is[[1]]), nrow = iter_mala, ncol = length(y))
  weights <- exp(as.numeric(unlist(mala.is[[2]])))
  
  
  asymp_covmat_fn <- function(chain, weights)
  {
  wts_mean <- mean(weights)
  num <- chain*weights
  sum_mat <- apply(num, 2, sum)
  is_est <- sum_mat / sum(weights)
  input_mat <- cbind(num, weights)  # input samples for mcse
  Sigma_mat <- mcse.multi(input_mat)$cov  # estimated covariance matrix of the tuple
  kappa_eta_mat <- cbind(diag(1/wts_mean, length(y)), -is_est/wts_mean) # derivative of kappa matrix
  asymp_covmat_is <- (kappa_eta_mat %*% Sigma_mat) %*% t(kappa_eta_mat) # IS asymptotic variance
  return(asymp_covmat_is)
  }
  
  #######  Quantiles for pi using Px
  ## MALA
  upper_quant_pi_lambda_pxm <- numeric(length = length(y))
  lower_quant_pi_lambda_pxm <- numeric(length = length(y))
  post_med_pi_lambda_pxm <- numeric(length = length(y))
  for(i in 1:length(y))
  {
    upper_quant_pi_lambda_pxm[i] <- quantile(output_single_run[[2]][,i], probs = 0.975)
    lower_quant_pi_lambda_pxm[i] <- quantile(output_single_run[[2]][,i], probs = 0.025)
    post_med_pi_lambda_pxm[i] <- quantile(output_single_run[[2]][,i], probs = 0.5)
  }
  
  
  
  
  
  pdf("plots/acf_tf_bark.pdf", height = 6, width = 12)
  ### Barker
  acf_isb <- acf(output_bark[,1], plot = FALSE)$acf
  acf_pxb <- acf(output_bark[,1], plot = FALSE)$acf
  
  plot(1:length(acf_isb), acf_isb, col = "blue", type = "l",
       xlab = "Lag", ylab = "Autocorrelation", ylim = c(-0.2, 1))
  lines(1:length(acf_pxb), acf_pxb, col = "red", type = "l")
  legend("bottomright", c("MYBarker", "PxBarker"), lty = 1,
         col = c("blue", "red"), cex = 0.75, bty = "n")
  
  for (i in 2:100) 
  {
    acf_isb <- acf(output_bark[,i], plot = FALSE)$acf
    acf_pxb <- acf(output_bark[,i], plot = FALSE)$acf
    lines(1:length(acf_isb), acf_isb, col = "blue", type = "l")
    lines(1:length(acf_pxb), acf_pxb, col = "red", type = "l")
  }
  dev.off()
  
  
  
  
  pdf("plots/acf_tf_bark.pdf", height = 6, width = 10)
  ### Barker
  acf_isb <- acf(mybark[[1]][,1], plot = FALSE)$acf
  acf_pxb <- acf(pxbark[[1]][,1], plot = FALSE)$acf
  
  plot(1:length(acf_isb), acf_isb, col = "blue", type = "l",
       xlab = "Lag", ylab = "Autocorrelation", ylim = c(-0.2, 1))
  lines(1:length(acf_pxb), acf_pxb, col = "red", type = "l")
  legend("bottomright", c("MYBarker", "PxBarker"), lty = 1,
         col = c("blue", "red"), cex = 0.75, bty = "n")
  
  for (i in 2:100) 
  {
    acf_isb <- acf(mybark[[1]][,i], plot = FALSE)$acf
    acf_pxb <- acf(pxbark[[1]][,i], plot = FALSE)$acf
    lines(1:length(acf_isb), acf_isb, col = "blue", type = "l")
    lines(1:length(acf_pxb), acf_pxb, col = "red", type = "l")
  }
  dev.off()
  
  
  
  
  
  # log_q_ratio_barker<-function(x,y,grad_x,grad_y)
  # {
  #   # x: current location (vector)
  #   # y: proposed location (vector)
  #   # grad_x: target log-posterior gradient at x (vector)
  #   # grad_y: target log-posterior gradient at y (vector)
  #   a1<-  c(-(x-y)*grad_y)
  #   a2<-  c(-grad_x*(y-x))
  #   return(sum(-(pmax(a1,0)+log1p(exp(-abs(a1))))+
  #                (pmax(a2,0)+log1p(exp(-abs(a2))))
  #   ))
  # }
  
