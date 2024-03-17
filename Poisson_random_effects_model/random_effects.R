####  Poisson random effects model
#### The function evaluated below of targets and the derivatives are omitting constant
#### of proportionality.

rm(list = ls())
set.seed(12345)
ni_s <- 5
I <- 50
c <- 10
sigma_eta <- 3
one_mat <- rep(1, I)
data <- matrix(0, nrow = I, ncol = ni_s)
identity_mat <- diag(1, I, I)
lambda <- 0.01
mu <- rnorm(1, 0, c)
eta_vec <- rnorm(I, mu, sigma_eta)
tol_nr <- 0.01
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
  term_exp_eta <- (ni_s*one_mat)*exp(eta)
  term_y <- apply(data, 1, sum)
  term_mu_sigma <- (mu/(sigma^2))*one_mat
  grad_vec <- - eta/(sigma^2) - term_exp_eta + 
    term_y + term_mu_sigma
  return(grad_vec)
}

true_hessian_eta <- function(sigma, eta)
{
  term_exp_diag <- 1/(sigma^2) + (ni_s*one_mat)*exp(eta)
  eta_hessian <- - diag(term_exp_diag)
}

proxfunc <- function(eta, mu, lambda, eta_initial, mu_initial, sigma)
{
  # For eta's
  grad_vec <- true_grad_vec(mu_initial, sigma, eta_initial) + (eta_initial - eta)/lambda
  eta_hessian <- true_hessian_eta(sigma, eta_initial) + 1/lambda
  while (sum(grad_vec^2) > tol_nr) 
    {
    eta_next <- eta_initial - solve(eta_hessian)%*%grad_vec
    grad_vec_next <- true_grad_vec(mu_initial, sigma, eta_next) + (eta_next - eta)/lambda
    eta_hessian_next <- true_hessian_eta(sigma, eta_next) + 1/lambda
    eta_initial <- eta_next
    grad_vec < grad_vec_next
  }
  mu_grad_vec <- mu_initial*I*(1/(sigma^2) + 1/(c^2)) - 
                           sum(abs(eta_initial))/(sigma^2) + (mu - mu_initial)/lambda
  mu_hessian <- I*(1/(sigma^2) + 1/(c^2)) + 1/lambda
  while (sum(mu_grad_vec^2) > tol_nr) 
    {
    mu_next <- mu_initial - mu_grad_vec/mu_hessian
    mu_grad_vec_next <- mu_next*I*(1/(sigma^2) + 1/(c^2)) - 
      sum(abs(eta_initial))/(sum(sigma^2)) + (mu - mu_next)/lambda
    mu_initial <- mu_next
    mu_grad_vec <- mu_
  }
  optima <- rbind(eta_initial, mu_initial)
  return(optima)
}

MY_env <- function(eta, mu, data, optima, lambda)
{
  term <- rbind(eta, mu)
  value <- true_target(eta, mu, data) + sum(term^2)/lambda
  return(-value)
}

eta_start <- eta_vec + rnorm(I, 0, 2)
mu_start <- 8
proxfunc(eta_vec, mu, lambda = lambda, eta_initial = eta_start, 
         mu_initial = mu_start, sigma = sigma_eta)

