####  Poisson random effects model
#### The function evaluated below of targets and the derivatives are omitting constant
#### of proportionality.

library(SimTools)
set.seed(12345)
ni_s <- 5
I <- 50
c <- 10
sigma_eta <- 3
one_mat <- rep(1, I)
data <- matrix(0, nrow = I, ncol = ni_s)
identity_mat <- diag(1, I, I)
lambda <- 0.001
mu <- rnorm(1, 0, c)
eta_vec <- rnorm(I, mu, sigma_eta)
tol_nr <- 1e-8
for (j in 1:ni_s)       # data generation
  {
  data[,j] <- rpois(I, exp(eta_vec))
}
data

true_target <- function(eta, mu, data)   # log of true target i.e. log(p(eta, mu|y))
{
  densval <- sum(eta^2)/(2*sigma_eta^2) - mu*sum(eta)/(sigma_eta^2) + 
                      ni_s*sum(exp(eta)) - sum(eta*apply(data, 1, sum)) +
                              (mu^2*I/2)*(1/(sigma_eta^2) + 1/(c^2))
  return(-densval)
}

true_grad_vec <- function(mu, sigma, eta)  # function evaluates gradient of log target
{
  term_exp_eta <- (ni_s)*exp(eta)
  term_y <- apply(data, 1, sum)
  term_mu_sigma <- (mu/(sigma^2))
  grad_vec_eta <- - eta/(sigma^2) - term_exp_eta + term_y + term_mu_sigma
  grad_mu <- - mu*(I/(sigma^2) + 1/(c^2)) + sum(eta)/(sigma^2)
  grad_value <- c(grad_vec_eta, grad_mu)
  return(grad_value)
}

true_hessian <- function(sigma, eta)  # function evaluates hessian of log target
{
  term_eta_diag <- 1/(sigma^2) + (ni_s)*exp(eta)
  # mu_term <- -(I/(sigma^2) + 1/(c^2))
  # mat_term_eta <- - diag(term_eta_diag, I, I)
  # vect <- -rep(1/(sigma^2), I)
  # vect_aug <- c(vect, mu_term)
  # mat_term_2 <- rbind(mat_term_eta, vect)
  # mat_final <- cbind(mat_term_2, vect_aug)

  return(-term_eta_diag)
}

proxfunc <- function(eta, mu, lambda, eta_initial, mu_initial, sigma)
{
  #  For starting values

  iter <- 0
  eta_next <- eta_initial
  mu_next <- mu_initial

  grad_vec <- true_grad_vec(mu_next, sigma, eta_next) - (c(eta_next, mu_next) - c(eta, mu))/lambda
  while (sqrt(sum(grad_vec^2)) > tol_nr) 
  {
     ## For eta's
    iter <- iter + 1

    # calculating hessian
    eta_hessian <- true_hessian(sigma, eta_next) - 1/lambda
    eta_grad <- grad_vec[1:I]    

    # defining the current state and the update 
    # current <- c(eta_next, mu_next)
    eta_next <- eta_next - eta_grad/eta_hessian


    mu.num <- (sum(eta_next)/sigma^2 + mu/lambda)
    mu.den <- I/(sigma^2) + 1/(c^2) + 1/lambda
    mu_next <- mu.num/mu.den
    # mu_next <- mu_initial - mu_grad/mu_hessian
    grad_vec <- true_grad_vec(mu_next, sigma, eta_next) - (c(eta_next, mu_next) - c(eta, mu))/lambda
  }
  #print(iter)
  optima <- c(eta_next, mu_next)
  return(optima)
}

MY_env <- function(eta, mu, data, optima, lambda)
{
  term <- rbind(eta, mu)
  value <- true_target(eta, mu, data) - sum(term^2)/(2*lambda)
  return(value)
}

log_gradpi <- function(eta, mu, lambda, eta_initial, mu_initial, sigma) #gradient of log target^lambda
{
  term_prox <- proxfunc(eta, mu, lambda, eta_initial, mu_initial, sigma)
  term <- c(eta, mu)
  ans <-  (term-term_prox)/lambda
  return(-ans)
}
eta_start <- eta_vec + rnorm(I, 0, 2)
mu_start <- 5
check <- proxfunc(eta_vec, mu, lambda = lambda, eta_initial = eta_vec, 
         mu_initial = mu, sigma = sigma_eta)

cbind(check, c(eta_vec,mu))

##### MYMALA samples function

mymala <- function(eta_start, mu_start, lambda, sigma, iter, delta, data)
{
  samp.mym <- matrix(0, nrow = iter, ncol = I+1)
  wts_is_est <- numeric(length = iter)
 
  samp_current <- c(eta_start, mu_start)
  samp.mym[1,] <- samp_current
  
  psi_val <- true_target(eta_start, mu_start, data)
  prox_val.curr <- proxfunc(eta_start, mu_start, lambda, eta_start, mu_start, sigma)
  psi_lambda_val <- true_target(prox_val.curr[1:I], prox_val.curr[I+1], data) -
    sum((prox_val.curr - c(eta_start, mu_start))^2)/(2*lambda)
  
  wts_is_est[1] <- psi_val - psi_lambda_val
 
   accept <- 0
  for (i in 2:iter) 
  {
    samp_next <- rnorm(length(samp_current), samp_current + 
                       (delta / 2)*log_gradpi(samp_current[1:I], samp_current[I+1],
                        lambda, eta_start, mu_start, sigma),  sqrt(delta))   # proposal step
    prox_val.next <- proxfunc(samp_next[1:I], samp_next[I+1],
                               lambda, eta_start, mu_start, sigma)
    
    targ_val.next <- true_target(prox_val.next[1:I], prox_val.next[I+1], data) -
                          sum((prox_val.next - c(eta_start, mu_start))^2)/(2*lambda)
    targ_val.curr <- true_target(prox_val.curr[1:I], prox_val.curr[I+1], data) -
                          sum((prox_val.curr - c(eta_start, mu_start))^2)/(2*lambda)
    
    q.next_to_curr <- sum(dnorm(samp_current, samp_next + 
                                  (delta / 2)*log_gradpi(samp_next[1:I], samp_next[I+1],
                                                         lambda, eta_start, mu_start, sigma),sqrt(delta), log = TRUE))
    q.curr_to_next <- sum(dnorm(samp_next, samp_current + 
                                  (delta / 2)*log_gradpi(samp_current[1:I], samp_current[I+1],
                                                         lambda, eta_start, mu_start, sigma),  sqrt(delta), log = TRUE))
    mh.ratio <- targ_val.next + q.next_to_curr - (targ_val.curr + q.curr_to_next)  # mh  ratio
    
    if(log(runif(1)) <= mh.ratio)
    {
      samp.mym[i,] <- samp_next
      prox_val.curr <- prox_val.next
      psi_val <- true_target(samp_next[1:I], samp_next[I+1], data)
      psi_lambda_val <- true_target(prox_val.next[1:I], prox_val.next[I+1], data) -
        sum((prox_val.next - samp_next)^2)/(2*lambda)
      wts_is_est[i] <-  psi_val - psi_lambda_val
      accept <- accept + 1
    }
    else
    {
      samp.mym[i,] <- samp_current
      psi_val <- true_target(samp_next[1:I], samp_next[I+1], data)
      psi_lambda_val <- true_target(prox_val.curr[1:I], prox_val.curr[I+1], data) -
        sum((prox_val.curr - samp_current)^2)/(2*lambda)
      wts_is_est[i] <- psi_val - psi_lambda_val
    }
    samp_current <- samp.mym[i,]
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(accept/iter)
  object <- list(samp.mym, wts_is_est)
  return(object)
}

##### PxMALA samples function

px.mala <- function(eta_start, mu_start, lambda, sigma, iter, delta, data)
{
  samp.pxm <- matrix(0, nrow = iter, ncol = I+1)
  samp_current <- c(eta_start, mu_start)
  samp.pxm[1,] <- samp_current
  accept <- 0
  
  for (i in 2:iter)
  {
    samp_next <- rnorm(length(samp_current), samp_current + 
                     (delta / 2)*log_gradpi(samp_current[1:I], samp_current[I+1],
                         lambda, eta_start, mu_start, sigma),  sqrt(delta))   # proposal step
    U_sampnext <- true_target(samp_next[1:I], samp_next[I+1], data)
    U_sampcurr <- true_target(samp_current[1:I], samp_current[I+1], data)
    q.next_to_curr <- sum(dnorm(samp_current, samp_next + 
                         (delta / 2)*log_gradpi(samp_next[1:I], samp_next[I+1],
                             lambda, eta_start, mu_start, sigma),sqrt(delta), log = TRUE))
    q.curr_to_next <- sum(dnorm(samp_next, samp_current + 
                        (delta / 2)*log_gradpi(samp_current[1:I], samp_current[I+1],
                             lambda, eta_start, mu_start, sigma),  sqrt(delta), log = TRUE))
    mh.ratio <- U_sampnext + q.next_to_curr - (U_sampcurr + q.curr_to_next)
    if(log(runif(1)) <= mh.ratio)
    {
      samp.pxm[i,] <- samp_next
      accept <- accept + 1
    }
    else
    {
      samp.pxm[i,] <- samp_current
    }
    samp_current <- samp.pxm[i,]
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(accept/iter)
  return(samp.pxm)
}


##  myhmc samples

myhmc <- function(eta_start, mu_start,lambda, sigma, iter, data, eps_hmc, L)
{
  samp.hmc <- matrix(0, nrow = iter, ncol = I+1)
  wts_is_est <- numeric(length = iter)
  samp <- c(eta_start, mu_start)
  samp.hmc[1,] <- samp
  psi_val <- true_target(eta_start, mu_start, data)
  proxval_curr <- proxfunc(eta_start, mu_start, lambda, eta_start, mu_start, sigma)
  psi_lambda_val <- true_target(proxval_curr[1:I], proxval_curr[I+1], data) +
    sum((proxval_curr - c(eta_start, mu_start))^2)/lambda
  wts_is_est[1] <- psi_lambda_val - psi_val
  mom_mat <- matrix(rnorm(iter*(I+1)), nrow = iter, ncol = I+1)
  accept <- 0
  
  for (i in 2:iter) 
  {
    p_prop <- mom_mat[i,]
    U_samp <- -log_gradpi(samp[1:I], samp[I+1],
                          lambda, eta_start, mu_start, sigma)
    p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
    q_current <- samp
    for (j in 1:L)
    {
      samp <- samp + eps_hmc*p_current   # full step for position
      U_samp <- -log_gradpi(samp[1:I], samp[I+1],
                            lambda, eta_start, mu_start, sigma)
      if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
    }
    p_current <- p_current - eps_hmc*U_samp/2
    p_current <- - p_current  # negation to make proposal symmetric
    
    #  proximal values
    proxval_prop <-  proxfunc(samp[1:I], samp[I+1], lambda, eta_start, mu_start, sigma)
    
    U_curr <- - (true_target(proxval_curr[1:I], proxval_curr[I+1], data) +
                   sum((proxval_curr - q_current)^2)/lambda)
    U_prop <- - (true_target(proxval_prop[1:I], proxval_prop[I+1], data) +
                   sum((proxval_curr - samp)^2)/lambda)
    K_curr <-  sum((p_prop^2)/2)
    K_prop <-  sum((p_current^2)/2)
    
    log_acc_prob = U_curr - U_prop + K_curr - K_prop
    
    if(log(runif(1)) <= log_acc_prob )
    {
      samp.hmc[i,] <- samp
      proxval_curr <- proxval_prop
      psi_val <- true_target(samp[1:I], samp[I+1], data)
      psi_lambda_val <- - (true_target(proxval_prop[1:I], proxval_prop[I+1], data) +
                             sum((proxval_prop - samp)^2)/lambda)
      wts_is_est[i] <- psi_lambda_val - psi_val
      accept <- accept + 1
    }
    else
    {
      samp.hmc[i,] <- q_current
      psi_val <- true_target(q_current[1:I], q_current[I+1], data)
      psi_lambda_val <- - (true_target(proxval_curr[1:I], proxval_curr[I+1], data) +
                             sum((proxval_curr - q_current)^2)/lambda)
      wts_is_est[i] <- psi_lambda_val - psi_val
      samp <- q_current
    }
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(acc_rate <- accept/iter)
  object <- list(samp.hmc, wts_is_est, acc_rate)
  return(object)
}

## pxhmc samples

pxhmc <- function(eta_start, mu_start,lambda, sigma, iter, data, eps_hmc, L)
{
  samp.hmc <- matrix(0, nrow = iter, ncol = I+1)
  samp <- c(eta_start, mu_start)
  samp.hmc[1,] <- samp
  mom_mat <- matrix(rnorm(iter*(I+1)), nrow = iter, ncol = I+1)
  accept <- 0
  
  for (i in 2:iter) 
  {
    p_prop <- mom_mat[i,]
    U_samp <- -log_gradpi(samp[1:I], samp[I+1],
                          lambda, eta_start, mu_start, sigma)
    p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
    q_current <- samp
    for (j in 1:L)
    {
      samp <- samp + eps_hmc*p_current   # full step for position
      U_samp <- -log_gradpi(samp[1:I], samp[I+1],
                            lambda, eta_start, mu_start, sigma)
      if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
    }
    p_current <- p_current - eps_hmc*U_samp/2
    p_current <- - p_current  # negation to make proposal symmetric
    
    U_curr <- - true_target(q_current[1:I], q_current[I+1], data)
    U_prop <- - true_target(samp[1:I], samp[I+1], data)
    K_curr <-  sum((p_prop^2)/2)
    K_prop <-  sum((p_current^2)/2)
    
    log_acc_prob = U_curr - U_prop + K_curr - K_prop
    
    if(log(runif(1)) <= log_acc_prob )
    {
      samp.hmc[i,] <- samp
      accept <- accept + 1
    }
    else
    {
      samp.hmc[i,] <- q_current
      samp <- q_current
    }
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(acc_rate <- accept/iter)
  object <- list(samp.hmc, acc_rate)
  return(object)
}

eta_start <- log(rowMeans(data))
mu_start <- mean(eta_start)
lambda <- 100
t <- mymala(eta_start, mu_start, lambda,sigma_eta, iter = 1e4, delta = 1.5, data)
u <- px.mala(eta_start, mu_start, lambda,sigma_eta, iter = 1e4, delta = 7e-7, data)
plot.ts(t[[2]])
library(SimTools)
plot.ts(t[[1]][,1:10])
plot.ts(t[[1]][,21:30])
plot.ts(t[[1]][,31:40])
plot.ts(t[[1]][,41:50])
plot.ts(t[[1]][,51])

acfplot(t[[1]][, 1:5])

mala_chain <- t[[1]]
weights <- t[[2]]
is_samp <- matrix(unlist(mala_chain), nrow = 1e4, ncol = I+1)
is_wts <- as.numeric(unlist(weights))
wts_mean <- mean(exp(is_wts))
num <- is_samp*exp(is_wts)
sum_mat <- apply(num, 2, sum)
is_est <- sum_mat / sum(exp(is_wts))
input_mat <- cbind(num, exp(is_wts))  # input samples for mcse
Sigma_mat <- mcse.multi(input_mat)$cov  # estimated covariance matrix of the tuple
kappa_eta_mat <- cbind(diag(1/wts_mean, I+1), is_est/wts_mean) # derivative of kappa matrix

asymp_covmat_is <- (kappa_eta_mat %*% Sigma_mat) %*% t(kappa_eta_mat) # IS asymptotic variance

asymp_covmat_pxm <- mcse.multi(u)$cov   # PxMALA asymptotic variance

rel_ess <- (det(asymp_covmat_pxm)/det(asymp_covmat_is))^(1/(I+1))

diag(asymp_covmat_pxm)/diag(asymp_covmat_is)


v <- myhmc(eta_start, mu_start,lambda, sigma_eta, iter = 1e3, data, eps_hmc=0.05, L=10)
w <- pxhmc(eta_start, mu_start,lambda, sigma_eta, iter = 1e3, data, eps_hmc=0.05, L=10)



