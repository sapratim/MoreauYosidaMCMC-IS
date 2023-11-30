### One dimensional Bayesian Lasso

beta_par <- 3  # true value of beta
ndata <- 100
noise_sd <- 3 
x <- seq(1, 10, length = 100)   # vector of regressors
y <- rnorm(length(x), x*beta_par, noise_sd)  # actual generated data
alpha_true <- 5  # penalty parameter
lamb_coeff <- 1  # coefficient of closeness
step_size <- 2
beta_start <- 1   
iterations <- 1e5  

log_post <- function(y,beta,alpha,sigma2)  # log value of the posterior distribution
 {
  term1 <- (sum((y-x*beta)^2))/(2*sigma2)
  term2 <- (alpha*abs(beta))/sqrt(sigma2)
  value <- term1 + term2
  return(-value)
}

log_target <- function(eta,beta,lambda,y,sigma2,alpha)  # log of MY envelope
{
  dens_val <- sum((y - x*eta)^2)/ (2*sigma2) + (alpha*abs(eta))/sqrt(sigma2) + 
                                  sum((eta - beta)^2)/(2*lambda) 
  return(-dens_val)
}

prox_func <- function(beta,lambda,y,sigma2,alpha)  # function returns proximal value
{
  denom <- lambda*(sum(x^2)) + sigma2
  mu <- (lambda*(sum(x*y)) + beta*sigma2) / denom 
  temp <- (alpha*lambda*sqrt(sigma2)) / denom
  vec <- c(0, mu + temp, mu - temp)
  fun.val <- numeric(length = length(vec))
  fun.val[1] <- log_target(vec[1],beta,lambda,y,sigma2,alpha)
  fun.val[2] <- log_target(vec[2],beta,lambda,y,sigma2,alpha)
  fun.val[3] <- log_target(vec[3],beta,lambda,y,sigma2,alpha)
  index <- which.min(fun.val)
  prox <- vec[index]
  return(prox)
}

log_gradpi <- function(beta,lambda,y,sigma2,alpha)  # gradient of log target(MY envelope)
{
  beta_prox <- prox_func(beta,lambda,y,sigma2,alpha)
  ans <-  (beta-beta_prox)/lambda
  return(-ans)
}

mymala_fn <- function(y, alpha, sigma2, lambda, iter, delta)
{
  samp.mym <- matrix(0, nrow = iter, ncol = 1)
  beta_current <- beta_start
  samp.mym[1,] <- beta_current
  wts_is_est <- numeric(length = iter)
  prox_beta_current <- prox_func(beta_current,lambda,y,sigma2,alpha)
  g_lambda_val <- log_target(prox_beta_current,beta_current,lambda,y,sigma2,alpha)
  g_val <- log_post(y,beta_current,alpha,sigma2)
  wts_is_est[1] <- g_lambda_val - g_val
  accept <- 0
  for (i in 2:iter) 
  {
    beta_next <- rnorm(length(beta_current), beta_current + 
                         (delta / 2)*log_gradpi(beta_current,lambda,y,sigma2,alpha), 
                       sqrt(delta))   # proposal step
    prox_val.next <- prox_func(beta_next, lambda, y, sigma2, alpha)
    prox_val.curr <- prox_func(beta_current, lambda, y, sigma2, alpha)
    targ_val.next <- log_target(prox_val.next,beta_next,lambda,y,sigma2,alpha)
    targ_val.curr <- log_target(prox_val.curr,beta_current,lambda,y,sigma2,alpha)
    q.next_to_curr <- sum(dnorm(beta_current, beta_next + 
                                  (delta / 2)*log_gradpi(beta_next,lambda,y,sigma2,alpha),
                                sqrt(delta), log = TRUE))
    q.curr_to_next <- sum(dnorm(beta_next, beta_current + 
                                  (delta / 2)*log_gradpi(beta_current,lambda,y,sigma2,alpha),
                                sqrt(delta), log = TRUE))
    mh.ratio <- targ_val.next + q.next_to_curr - (targ_val.curr + q.curr_to_next)  # mh  ratio
    # print(mh.ratio)
    if(log(runif(1)) <= mh.ratio)
    {
      samp.mym[i,] <- beta_next
      prox_beta_next <- prox_func(beta_next,lambda,y,sigma2,alpha)
      g_val <- log_post(y,beta_next,alpha,sigma2)
      g_lambda_val <- log_target(prox_beta_next,beta_next,lambda,y,sigma2,alpha)
      wts_is_est[i] <- g_lambda_val - g_val
      accept <- accept + 1
    }
    else
    {
      samp.mym[i,] <- beta_current
      prox_beta_current <- prox_func(beta_current,lambda,y,sigma2,alpha)
      g_val <- log_post(y,beta_current,alpha,sigma2)
      g_lambda_val <- log_target(prox_beta_current,beta_current,lambda,y,sigma2,alpha)
      wts_is_est[i] <- g_lambda_val - g_val
    }
    beta_current <- samp.mym[i,]
    if(i %% 10000 == 0){
      print(i)
    }
  }
  print(accept/iter)
  object <- list(samp.mym, wts_is_est)
  return(object)
}

result <- mymala_fn(y=y,alpha = alpha_true,sigma2 = noise_sd^2,lambda = lamb_coeff, 
               iter = iterations, delta = step_size)

# Mode matching
values <- seq(-35,40, length = 1e4)
prox_vals <- numeric(length = 1e4)
den_val <- numeric(length = 1e4)
for (i in 1:length(prox_vals)) 
  {
  prox_vals[i] = prox_func(values[i],lamb_coeff,y,noise_sd^2,alpha_true)
  den_val[i] <- log_target(prox_vals[i],values[i],lamb_coeff,y,noise_sd^2,alpha_true)
}
plot(values, den_val, type = "l")
abline(v=prox_func(beta=3,lamb_coeff,y,noise_sd^2,alpha_true), col = "red")
abline(v=3, col = "blue")
