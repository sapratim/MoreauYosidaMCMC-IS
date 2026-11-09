# Copies of the original functions, with proposal records added.
# Sampling equations, caching, random draws and acceptance tests are unchanged.
# The first proposal/ADMM/accepted entries are NA: state 1 is initialization.
# Relevant functions copied from ../TF_functions.R.
# Only the HMC code and its dependencies are needed for this check.
library(glmgen)
library(Matrix)

# generating data
set.seed(12345)
alpha_hat <- 5   # obtained from the first dataset
sigma2_hat <- 9  # obtained from the first dataset
x <- seq(1,100,len=100)
f <- Vectorize(function(x){if(x<=35){x} else if(x<=70){70-x} else{0.5*x-35}})
fx_linear <- f(x)
y <- fx_linear + rnorm(length(x), sd = 3)


# setting parameters for the proximal mappings
tol <- 1e-5
max_iter <- 200L
trend_control <- trendfilter.control.list(obj_tol = tol, max_iter = max_iter)

# function calculates the D matrix of the penalty term
getD <- function(k, n, x=NULL){
  if(is.null(x)){
    x <- 1:n
  }
  diags <- list(rep(-1,n),rep(1,n))
  D <- Matrix::bandSparse(n-1,n,k=c(0,1),diag=diags,symm=F)
  if(k>=1){
    for(i in 1:k){
      leftD <- Matrix::bandSparse(n-i-1,n-i,k=c(0,1),diag=diags,symm=F)
      xdiag <- Matrix::Diagonal(n-i,i/diff(x,lag=i))
      D <- leftD %*% xdiag %*% D
    }
  }
  return(D)
}

# log target of pi-lambda
log_pilambda <- function(eta,beta,lambda,y,sigma2,alpha)
{
  dens_val <- alpha*(sum(abs(D_mat%*%eta))) + sum((beta-eta)^2)/(2*lambda) + 
    sum((y-eta)^2)/(2*sigma2)
  return(-dens_val)
}

log_pi <- function(beta,y,sigma2,alpha)
{
  dens_val <- alpha*(sum(abs(D_mat%*%beta))) + sum((y-beta)^2)/(2*sigma2)
  return(-dens_val)
}

# Same proximal mapping; also return glmgen's reported ADMM count.
prox_record <- function(beta, lambda, alpha, sigma2, k, grid)
{
  betaval <- (beta*sigma2 + lambda*y)/(sigma2 + lambda)
  lambdaval <- alpha / ((lambda + sigma2)/(lambda*sigma2))
  fit <- trendfilter(grid, betaval, k = k, lambda = lambdaval,
                     control = trend_control)
  list(beta = as.vector(fit$beta), admm_iterations = fit$iter)
}

##  myhmc samples

myhmc_record <- function(y, alpha, sigma2, k, grid, iter, eps_hmc, L, start)
{
  nvar <- length(y)
  samp.hmc <- matrix(0, nrow = iter, ncol = nvar)
  proposal_states <- matrix(NA_real_, nrow = iter, ncol = nvar)
  proposal_admm <- rep(NA_integer_, iter)
  accepted <- rep(NA, iter)
  lambda <- lamb_coeff
  wts_is_est <- numeric(iter)
  L_val <- L
  
  ## initial state
  beta <- start
  samp.hmc[1, ] <- beta
  
  ## cache prox and gradient at current state
  prox_fit <- prox_record(beta, lambda, alpha, sigma2, k, grid)
  proxval_curr <- prox_fit$beta
  grad_curr <- -(beta - proxval_curr)/lambda
  
  ## current MY potential
  U_curr <- -log_pilambda(
    proxval_curr, beta, lambda, y, sigma2, alpha
  )
  
  ## initial weight
  g_val <- -log_pi(beta, y, sigma2, alpha)
  wts_is_est[1] <- U_curr - g_val
  
  #mom_mat <- matrix(rnorm(iter*nvar), nrow = iter, ncol = nvar)
  accept <- 0
  
  for (i in 2:iter)
  {
    p_prop <- rnorm(nvar) #mom_mat[i, ]
    q_current <- beta
    
    ## use cached current gradient
    U_beta <- -grad_curr
    p_current <- p_prop - eps_hmc * U_beta / 2
    L_run <- ifelse(runif(1) <= 0.05, 1, L_val)
    
    for (j in 1:L_run)
    {
      beta <- beta + eps_hmc * p_current
      
      ## one prox evaluation at each new position
      prox_fit <- prox_record(beta, lambda, alpha, sigma2, k, grid)
      prox_beta <- prox_fit$beta
      grad_beta <- -(beta - prox_beta)/lambda
      U_beta <- -grad_beta
      
      if (j != L_run)
        p_current <- p_current - eps_hmc * U_beta
    }
    
    p_current <- p_current - eps_hmc * U_beta / 2
    p_current <- -p_current

    # Record the endpoint before the accept/reject decision.
    proposal_states[i, ] <- beta
    proposal_admm[i] <- prox_fit$admm_iterations
    
    ## final proposal prox is already prox_beta
    proxval_prop <- prox_beta
    grad_prop <- grad_beta
    
    U_prop <- -log_pilambda(
      proxval_prop, beta, lambda, y, sigma2, alpha
    )
    
    K_curr <- sum(p_prop^2)/2
    K_prop <- sum(p_current^2)/2
    
    log_acc_prob <- U_curr - U_prop + K_curr - K_prop
    
    if (log(runif(1)) <= log_acc_prob)
    {
      accepted[i] <- TRUE
      samp.hmc[i, ] <- beta
      
      ## cache accepted state
      proxval_curr <- proxval_prop
      grad_curr <- grad_prop
      U_curr <- U_prop
      
      ## weight
      g_val <- -log_pi(beta, y, sigma2, alpha)
      wts_is_est[i] <- U_curr - g_val
      
      accept <- accept + 1
    }
    else
    {
      accepted[i] <- FALSE
      beta <- q_current
      samp.hmc[i, ] <- beta
      
      ## all cached quantities remain unchanged
      wts_is_est[i] <- U_curr - g_val
    }
    
    if (i %% (iter/10) == 0) {
      print(cat(i, accept/iter))
    }
  }
  
  acc_rate <- accept/iter
  print(acc_rate)
  
  list(samples = samp.hmc, log_weights = wts_is_est, acceptance = acc_rate,
       proposals = proposal_states, admm_iterations = proposal_admm, accepted = accepted)
}

## pxhmc samples
pxhmc_record <- function(y, alpha, sigma2, k, grid, iter, eps_hmc, L, start)
{
  nvar <- length(y)
  samp.hmc <- matrix(0, nrow = iter, ncol = nvar)
  proposal_states <- matrix(NA_real_, nrow = iter, ncol = nvar)
  proposal_admm <- rep(NA_integer_, iter)
  accepted <- rep(NA, iter)
  lambda <- lamb_coeff
  L_val <- L
  
  beta <- start
  samp.hmc[1, ] <- beta
  
  ## cache current gradient
  prox_fit <- prox_record(beta, lambda, alpha, sigma2, k, grid)
  proxval_curr <- prox_fit$beta
  grad_curr <- -(beta - proxval_curr)/lambda
  
  ## cache current target value
  U_curr <- -log_pi(beta, y, sigma2, alpha)
  
  #mom_mat <- matrix(rnorm(iter*nvar), nrow = iter, ncol = nvar)
  accept <- 0
  
  for (i in 2:iter)
  {
    p_prop <- rnorm(nvar) #mom_mat[i, ]
    q_current <- beta
    
    ## cached gradient
    U_beta <- -grad_curr
    p_current <- p_prop - eps_hmc * U_beta / 2
    L_run <- ifelse(runif(1) <= 0.05, 1, L_val)
    
    for (j in 1:L_run)
    {
      beta <- beta + eps_hmc * p_current
      
      prox_fit <- prox_record(beta, lambda, alpha, sigma2, k, grid)
      prox_beta <- prox_fit$beta
      grad_beta <- -(beta - prox_beta)/lambda
      U_beta <- -grad_beta
      
      if (j != L_run)
        p_current <- p_current - eps_hmc * U_beta
    }
    
    p_current <- p_current - eps_hmc * U_beta / 2
    p_current <- -p_current

    # Record the endpoint before the accept/reject decision.
    proposal_states[i, ] <- beta
    proposal_admm[i] <- prox_fit$admm_iterations
    
    ## only target at proposal needs evaluating
    U_prop <- -log_pi(beta, y, sigma2, alpha)
    
    K_curr <- sum(p_prop^2)/2
    K_prop <- sum(p_current^2)/2
    
    log_acc_prob <- U_curr - U_prop + K_curr - K_prop
    
    if (log(runif(1)) <= log_acc_prob)
    {
      accepted[i] <- TRUE
      samp.hmc[i, ] <- beta
      
      ## cache everything at accepted proposal
      proxval_curr <- prox_beta
      grad_curr <- grad_beta
      U_curr <- U_prop
      
      accept <- accept + 1
    }
    else
    {
      accepted[i] <- FALSE
      beta <- q_current
      samp.hmc[i, ] <- beta
      
      ## proxval_curr, grad_curr and U_curr remain unchanged
    }
    
    if (i %% (iter/10) == 0) {
      print(cat(i, accept/iter))
    }
  }
  
  acc_rate <- accept/iter
  print(acc_rate)
  
  list(samples = samp.hmc, acceptance = acc_rate,
       proposals = proposal_states, admm_iterations = proposal_admm, accepted = accepted)
}
