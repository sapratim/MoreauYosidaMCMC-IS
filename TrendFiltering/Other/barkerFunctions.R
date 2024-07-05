### Barker's proposal functions

bark.prop <- function(beta,lambda,y,sigma2,alpha,k,grid,delta,covmat)
{
  aux_var <- rnorm(length(beta), 0, 1)
  y <- (delta*covmat)%*%aux_var
  denom_prod <- y*log_gradpi(beta,lambda,y,sigma2,alpha,k,grid)
  prob <- 1 / (1 + exp(- sum(denom_prod)))
  ifelse(runif(1) <= prob, prop <- beta + y, prop <- beta - y)
  return(prop)
}

log_q_ratio_barker<-function(x,y,grad_x,grad_y)
{
  # x: current location (vector)
  # y: proposed location (vector)
  # grad_x: target log-posterior gradient at x (vector)
  # grad_y: target log-posterior gradient at y (vector)
  a1<-  c(-(x-y)*grad_y)
  a2<-  c(-grad_x*(y-x))
  return(sum(-(pmax(a1,0)+log1p(exp(-abs(a1))))+
               (pmax(a2,0)+log1p(exp(-abs(a2))))
  ))
}


## IS from Barker

mybarker <- function(y, alpha, sigma2, k, grid, iter, delta, covmat)
{
  samp.bark <- matrix(0, nrow = iter, ncol = length(y))
  lambda <- lamb_coeff
  wts_is_est <- numeric(length = iter)
  beta_current <- markov_chain[1e5,]
  samp.bark[1,] <- beta_current
  g_val <- alpha*sum(abs(D_mat%*%beta_current)) + sum((y - beta_current)^2)/(2*sigma2)
  prox_val.curr <- prox_func(beta_current, lambda, alpha, sigma2, k, grid)
  g_lambda_val <- prox_arg(prox_val.curr, beta_current, lambda=lambda, y, sigma2, alpha)
  wts_is_est[1] <- g_lambda_val - g_val
  U <- sqrtm(covmat)
  mat.inv <- solve(covmat)
  accept <- 0
  for (i in 2:iter) 
  {
    beta_next <- bark.prop(beta_current,lambda,y,sigma2,alpha,k,grid,delta,covmat) # proposal step
    prox_val.next <- prox_func(beta_next, lambda, alpha, sigma2, k, grid)
    targ_val.next <- log_target(prox_val.next,beta_next,lambda,y,sigma2,alpha)
    targ_val.curr <- log_target(prox_val.curr,beta_current,lambda,y,sigma2,alpha)
    grad_beta_curr <- log_gradpi(beta_current, lambda, y, sigma2, alpha, k, grid)
    grad_beta_next <- log_gradpi(beta_next, lambda, y, sigma2, alpha, k, grid)
    
    mh.ratio <- targ_val.next - targ_val.curr + 
      log_q_ratio_barker(beta_current,beta_next,grad_beta_curr,grad_beta_next)
    
    if(log(runif(1)) <= mh.ratio)
    {
      samp.bark[i,] <- beta_next
      prox_val.curr <- prox_val.next
      g_val <- alpha*sum(abs(D_mat%*%beta_next)) + sum((y - beta_next)^2)/(2*sigma2)
      g_lambda_val <- prox_arg(prox_val.next, beta_next, lambda=lambda, y, sigma2, alpha)
      wts_is_est[i] <- g_lambda_val - g_val
      accept <- accept + 1
    }
    else
    {
      samp.bark[i,] <- beta_current
      g_val <- alpha*sum(abs(D_mat%*%beta_current)) + sum((y - beta_current)^2)/(2*sigma2)
      g_lambda_val <- prox_arg(prox_val.curr, beta_current, lambda=lambda, y, sigma2, alpha)
      wts_is_est[i] <- g_lambda_val - g_val
    }
    beta_current <- samp.bark[i,]
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(acc_rate <- accept/iter)
  object <- list(samp.bark, wts_is_est, acc_rate)
  return(object)
}

## PxBarker samples

px.barker <- function(y, alpha, sigma2, k, grid, iter, delta, covmat)
{
  samp.bark <- matrix(0, nrow = iter, ncol = length(y))
  lambda <- lamb_coeff
  beta_current <- markov_chain[1e5,]
  samp.bark[1,] <- beta_current
  mat.inv <- solve(covmat)
  accept <- 0
  for (i in 2:iter)
  {
    beta_next <- bark.prop(beta_current,lambda,y,sigma2,alpha,k,grid,delta,covmat)
    U_betanext <- - (sum((y - beta_next)^2)/(2*sigma2) + alpha*(sum(abs(D_mat%*%beta_next))))
    U_betacurr <- - (sum((y - beta_current)^2)/(2*sigma2) + alpha*(sum(abs(D_mat%*%beta_current))))
    grad_beta_curr <- log_gradpi(beta_current, lambda, y, sigma2, alpha, k, grid)
    grad_beta_next <- log_gradpi(beta_next, lambda, y, sigma2, alpha, k, grid)
    
    mh.ratio <- U_betanext - U_betacurr + 
      log_q_ratio_barker(beta_current,beta_next,grad_beta_curr,grad_beta_next)
    if(log(runif(1)) <= mh.ratio)
    {
      samp.bark[i,] <- beta_next
      accept <- accept + 1
    }
    else
    {
      samp.bark[i,] <- beta_current
    }
    beta_current <- samp.bark[i,]
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(acc_rate <- accept/iter)
  object <- list(samp.bark, acc_rate)
  return(object)
}