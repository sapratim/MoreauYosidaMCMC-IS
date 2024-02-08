  # R Code for Proximal of the Nuclear-Norm problem:
    # softthreshold(seq(-10, 10, length.out=10), 1) 

rm(list = ls())
library(mcmcse)
library(coda)
library(Matrix)
library(expm)
library(foreach)
library(doParallel)
library(ks)

set.seed(11111)
n <- 8
a <- 2
mat <- matrix(0, nrow = n, ncol = n)
vec.mat <- rep(c(1, 0), each = a, times = n/(2*a))
checker <- matrix(0, nrow = n, ncol = n)
for(j in seq(1, n/2, by = a))
{
  for(k in 1:a)
  {
    checker[j+k-1, ] <- vec.mat
  }
  vec.mat <- rev(vec.mat)
}
for(j in seq((n/2+1), n, by = a))
{
  for(k in 1:a)
  {
    checker[j+k-1, ] <- (vec.mat == 0)*(vec.mat + .25) + (vec.mat == 1)*(vec.mat - .25)
  }
  vec.mat <- rev(vec.mat)
}

noise <- matrix(rnorm(n^2, 0, sqrt(0.01)), nrow = n, ncol = n)
image_mat <- checker + noise

x <- vec(checker)
y <- vec(image_mat)
######### Display the checkerboard matrix as an image

# image(checker, col = gray.colors(4, start = 0, end = 1), axes = FALSE)
# image(image_mat, col = gray.colors(4, start = 0, end = 1), axes = FALSE)

##############------Functions------##################################

nucl_norm <- function(vect)    ## vector input
{
  A <- matrix(vect, nrow = n, ncol = n)
  norm_val <- sum(svd(A)$d)
  return(norm_val)
}

log_target <- function(eta,x,lambda,y,sigma2,alpha)   # log target function of pi^lambda
{
  n_norm <- nucl_norm(eta)
  dens_val <- alpha*n_norm + sum((eta-x)^2)/(2*lambda) + 
                                            sum((y-eta)^2)/(2*sigma2)
  return(-dens_val)
}

#########  Soft threshold function

softthreshold <- function(u, lambda) {       ####  u is a vector
  return(sign(u)*sapply(u, FUN=function(x) {max(abs(u)-lambda,0)}))
}

######### Proximity mapping

prox_func <- function(x,lambda,y,sigma2,alpha) {   #### input x and y as a vector
  num <- lambda*y + sigma2*x
  denom <- lambda + sigma2
  mat <- matrix(num/denom, nrow = n, ncol = n)
  svdsol <- svd(mat)
  s <- softthreshold(svdsol$d, (alpha*sigma2*lambda)/denom)
  output <- svdsol$u %*% (s*t(svdsol$v))  # Multiply each row of V’ by singular values
  return(vec(output))
}

log_gradpi <- function(x,lambda,y,sigma2,alpha)  # gradient of log target
{
  x_prox <- prox_func(x,lambda,y,sigma2,alpha)
  ans <-  (x-x_prox)/lambda
  return(-ans)
}

# dmvnorm_fn <- function(point, mu, mat, delta)
# {
#   diff <- point - mu
#   exp_term <- (t(diff) %*% mat) %*% diff
#   den_value <- -(exp_term/(2*delta))
# }

#### MYMALA sampling for covariance matrix estimation

mymala_cov_fn <- function(y, alpha, sigma2, iter, delta)  #pre-conditioned mala
{
  samp.mym <- matrix(0, nrow = iter, ncol = length(y))
  lambda <- lamb_coeff
  samp_current <- y
  prox_val.curr <- prox_func(y,lambda,y,sigma2,alpha)
  samp.mym[1,] <- samp_current
  accept <- 0
  for (i in 2:iter) 
  {
    samp_next <- rnorm(length(samp_current), samp_current + 
                         (delta / 2)*log_gradpi(samp_current,lambda,y,sigma2,alpha), 
                       sqrt(delta))   # proposal step
    prox_val.next <- prox_func(samp_next,lambda,y,sigma2,alpha)
    targ_val.next <- log_target(prox_val.next,samp_next,lambda,y,sigma2,alpha)
    targ_val.curr <- log_target(prox_val.curr,samp_current,lambda,y,sigma2,alpha)
    q.next_to_curr <- sum(dnorm(samp_current, samp_next + 
                                  (delta / 2)*log_gradpi(samp_next,lambda,y,sigma2,alpha),
                                sqrt(delta), log = TRUE))
    q.curr_to_next <- sum(dnorm(samp_next, samp_current + 
                                  (delta / 2)*log_gradpi(samp_current,lambda,y,sigma2,alpha),
                                sqrt(delta), log = TRUE))
    mh.ratio <- targ_val.next + q.next_to_curr - (targ_val.curr + q.curr_to_next)  # mh  ratio
    # print(mh.ratio)
    if(log(runif(1)) <= mh.ratio)
    {
      samp.mym[i,] <- samp_next
      prox_val.curr <- prox_val.next
      accept <- accept + 1
    }
    else
    {
      samp.mym[i,] <- samp_current
    }
    samp_current <- samp.mym[i,]
    if(i %% 1000 == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(accept/iter)
  object <- samp.mym
  return(object)
}

##### MYMALA samples function

mymala <- function(y, alpha, sigma2, iter, delta)
{
  samp.mym <- matrix(0, nrow = iter, ncol = length(y))
  lambda <- lamb_coeff
  wts_is_est <- numeric(length = iter)
  samp_current <- y
  samp.mym[1,] <- samp_current
  g_val <- alpha*nucl_norm(samp_current) + sum((y - samp_current)^2)/(2*sigma2)
  prox_val.curr <- prox_func(samp_current, lambda, y, sigma2, alpha)
  g_lambda_val <- - log_target(prox_val.curr,samp_current,lambda,y,sigma2,alpha)
  wts_is_est[1] <- g_lambda_val - g_val
 # U <- sqrtm(covmat)
#  mat.inv <- solve(covmat)
  accept <- 0
  for (i in 2:iter) 
  {
    samp_next <- rnorm(length(samp_current), samp_current + 
                         (delta / 2)*log_gradpi(samp_current,lambda,y,sigma2,alpha), 
                       sqrt(delta))   # proposal step
    prox_val.next <- prox_func(samp_next,lambda,y,sigma2,alpha)
    targ_val.next <- log_target(prox_val.next,samp_next,lambda,y,sigma2,alpha)
    targ_val.curr <- log_target(prox_val.curr,samp_current,lambda,y,sigma2,alpha)
    q.next_to_curr <- sum(dnorm(samp_current, samp_next + 
                                  (delta / 2)*log_gradpi(samp_next,lambda,y,sigma2,alpha),
                                sqrt(delta), log = TRUE))
    q.curr_to_next <- sum(dnorm(samp_next, samp_current + 
                                  (delta / 2)*log_gradpi(samp_current,lambda,y,sigma2,alpha),
                                sqrt(delta), log = TRUE))
    mh.ratio <- targ_val.next + q.next_to_curr - (targ_val.curr + q.curr_to_next)  # mh  ratio
    
    if(log(runif(1)) <= mh.ratio)
    {
      samp.mym[i,] <- samp_next
      prox_val.curr <- prox_val.next
      g_val <- alpha*nucl_norm(samp_next) + sum((y - samp_next)^2)/(2*sigma2)
      g_lambda_val <- - log_target(prox_val.next,samp_next,lambda,y,sigma2,alpha)
      wts_is_est[i] <- g_lambda_val - g_val
      accept <- accept + 1
    }
    else
    {
      samp.mym[i,] <- samp_current
      g_val <- alpha*nucl_norm(samp_current) + sum((y - samp_current)^2)/(2*sigma2)
      g_lambda_val <- - log_target(prox_val.curr,samp_current,lambda,y,sigma2,alpha)
      wts_is_est[i] <- g_lambda_val - g_val
    }
    samp_current <- samp.mym[i,]
    if(i %% 1000 == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(accept/iter)
  object <- list(samp.mym, wts_is_est)
  return(object)
}

##### PxMALA samples function

px.mala <- function(y, alpha, sigma2, iter, delta)
{
  samp.pxm <- matrix(0, nrow = iter, ncol = length(y))
  lambda <- lamb_coeff
  samp_current <- y
  samp.pxm[1,] <- samp_current
  accept <- 0
  # U <- sqrtm(covmat)
  # mat.inv <- solve(covmat)
  for (i in 2:iter)
  {
    samp_next <- rnorm(length(samp_current), samp_current + 
                         (delta / 2)*log_gradpi(samp_current,lambda,y,sigma2,alpha), 
                       sqrt(delta))   # proposal step
    U_sampnext <- - (sum((y - samp_next)^2)/(2*sigma2) + alpha*nucl_norm(samp_next))
    U_sampcurr <- - (sum((y - samp_current)^2)/(2*sigma2) + alpha*nucl_norm(samp_current))
    q.next_to_curr <- sum(dnorm(samp_current, samp_next + 
                                  (delta / 2)*log_gradpi(samp_next,lambda,y,sigma2,alpha),
                                sqrt(delta), log = TRUE))
    q.curr_to_next <- sum(dnorm(samp_next, samp_current + 
                                  (delta / 2)*log_gradpi(samp_current,lambda,y,sigma2,alpha),
                                sqrt(delta), log = TRUE))
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
    if(i %% 1000 == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(accept/iter)
  return(samp.pxm)
}


##  myhmc samples

myhmc <- function(y, alpha, sigma2, iter, eps_hmc, L)
{
  samp.hmc <- matrix(0, nrow = iter, ncol = length(y))
  lambda <- lamb_coeff
  wts_is_est <- numeric(length = iter)
  samp <- y
  samp.hmc[1,] <- y
  g_val <- alpha*nucl_norm(samp) + sum((y - samp)^2)/(2*sigma2)
  proxval_curr <- prox_func(samp, lambda, y, sigma2, alpha)
  g_lambda_val <- - log_target(proxval_curr,samp,lambda,y,sigma2,alpha)
  wts_is_est[1] <- g_lambda_val - g_val
  mom_mat <- matrix(rnorm(iter*length(y)), nrow = iter, ncol = length(y))
  accept <- 0
  
  for (i in 2:iter) 
  {
    p_prop <- mom_mat[i,]
    U_samp <- -log_gradpi(samp, lambda,y,sigma2,alpha)
    p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
    q_current <- samp
    for (j in 1:L)
    {
      samp <- samp + eps_hmc*p_current   # full step for position
      U_samp <- -log_gradpi(samp, lambda,y,sigma2,alpha)
      if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
    }
    p_current <- p_current - eps_hmc*U_samp/2
    p_current <- - p_current  # negation to make proposal symmetric
    
    #  proximal values
    proxval_prop <- prox_func(samp,lambda,y,sigma2,alpha)
    
    U_curr <- - log_target(proxval_curr,q_current,lambda,y,sigma2,alpha)
    U_prop <- - log_target(proxval_prop,samp,lambda,y,sigma2,alpha)
    K_curr <-  sum((p_prop^2)/2)
    K_prop <-  sum((p_current^2)/2)
    
    log_acc_prob = U_curr - U_prop + K_curr - K_prop
    
    if(log(runif(1)) <= log_acc_prob )
    {
      samp.hmc[i,] <- samp
      proxval_curr <- proxval_prop
      g_val <- alpha*nucl_norm(samp) + sum((y - samp)^2)/(2*sigma2)
      g_lambda_val <- - log_target(proxval_prop,samp,lambda,y,sigma2,alpha)
      wts_is_est[i] <- g_lambda_val - g_val
      accept <- accept + 1
    }
    else
    {
      samp.hmc[i,] <- q_current
      g_val <- alpha*nucl_norm(q_current) + sum((y - q_current)^2)/(2*sigma2)
      g_lambda_val <- - log_target(proxval_curr,samp,lambda,y,sigma2,alpha)
      wts_is_est[i] <- g_lambda_val - g_val
      samp <- q_current
    }
    if(i %% 1000 == 0){
      j <- accept/iter
      print(cat(i, j))}
  } 
  print(acc_rate <- accept/iter)
  object <- list(samp.hmc, wts_is_est, acc_rate)
  return(object)
}

## pxhmc samples

pxhmc <- function(y, alpha, sigma2, iter, eps_hmc, L)
{
  samp.hmc <- matrix(0, nrow = iter, ncol = length(y))
  lambda <- lamb_coeff
  samp <- y
  samp.hmc[1,] <- samp
  mom_mat <- matrix(rnorm(iter*length(y)), nrow = iter, ncol = length(y))
  accept <- 0
  
  for (i in 2:iter) 
  {
    p_prop <- mom_mat[i,]
    U_samp <- -log_gradpi(samp, lambda,y,sigma2,alpha)
    p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
    q_current <- samp
    for (j in 1:L)
    {
      samp <- samp + eps_hmc*p_current   # full step for position
      U_samp <- -log_gradpi(samp, lambda,y,sigma2,alpha)
      if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
    }
    p_current <- p_current - eps_hmc*U_samp/2
    p_current <- - p_current  # negation to make proposal symmetric
    
    U_curr <- sum((y - q_current)^2)/(2*sigma2) + alpha*nucl_norm(q_current)
    U_prop <- sum((y - samp)^2)/(2*sigma2) + alpha*nucl_norm(samp)
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
    if(i %% 1000 == 0){
      j <- accept/iter
      print(cat(i, j))}
  } 
  print(acc_rate <- accept/iter)
  object <- list(samp.hmc, acc_rate)
  return(object)
}


##########################################################################
##################  Program run ##########################################
##########################################################################


iter <- 1e4
lamb_coeff <- 0.00001
sigma2_hat <- 0.01
alpha_hat <- 1.15/sigma2_hat
step <- 0.0000001
step_pxm <- 0.0000005

result_is <- mymala(y=y, alpha = alpha_hat, sigma2 = sigma2_hat, iter = iter, delta = step)
result_pxm <- px.mala(y=y, alpha = alpha_hat, sigma2 = sigma2_hat, iter = iter, delta = step_pxm)

result_ishmc <- myhmc(y=y, alpha = alpha_hat, sigma2 = sigma2_hat, 
                               iter = iter, eps_hmc = 0.0001, L=10)
result_pxhmc <- pxhmc(y=y, alpha = alpha_hat, sigma2 = sigma2_hat, 
                      iter = iter, eps_hmc = 0.0001, L=10)

chain <- result_is[[1]]
weights <- result_is[[2]]
hist(weights)


############# Posterior mean

weight_mat <- matrix(0, nrow = 1e4, ncol = length(y))
for (i in 1:1e4) {
  weight_mat[i,] <- chain[i,]*exp(weights[i])
}
num_sum <- apply(weight_mat, 2, sum)
weights_sum <- sum(exp(weights))
post_mean <- num_sum/weights_sum


############ Image visualisation
chain1 <- result_pxm
image(matrix(post_mean, nrow = n, ncol = n), col = gray.colors(4, start = 0, end = 1), axes = FALSE)
library(SimTools)
traceplot(chain[, c(60:64)])
