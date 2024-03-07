####  Poisson random effects model

rm(list = ls())
set.seed(12345)
ni_s <- 5
I <- 50
c <- 10
sigma_eta <- 3
one_mat <- rep(1, I)
data <- matrix(0, nrow = I, ncol = ni_s)
identity_mat <- diag(1, I, I)
mu <- rnorm(1, 0, c)
for (j in 1:ni_s)       # data generation
  {
  eta_vec <- rnorm(I, mu, sigma_eta)
  data[,j] <- rpois(I, exp(eta_vec))
}
data

true_target <- function(eta, mu, data)
{
  densval <- sum(eta^2)/(2*sigma_eta^2) - mu*sum(eta)/(sigma_eta^2) + 
                      ni_s*sum(exp(eta)) - sum(eta*apply(data, 1, sum)) +
                              (mu^2*I/2)*(1/(sigma_eta^2) + 1/(c^2))
  return(-densval)
}

proxfunc <- function(lambda, eta_initial, mu_initial)
{
  # For eta's
  term_exp_eta <- (ni_s*t(one_mat))%*%exp(eta_initial)
  term_y <- apply(data, 1, sum)
  term_mu_sigma <- (mu/(sigma_eta^2))*one_mat
  grad_vec <- - eta_initial/(sigma_eta^2) - term_exp_eta + 
                term_y + term_mu_sigma
  eta_hessian <- - (1/(sigma_eta^2))%*%identity_mat - diag(term_exp_eta)
  while (sum(grad_vec^2) > 0.0001) 
    {
    eta_next <- eta_initial - solve(eta_hessian)%*%grad_vec
    term_exp_eta <- (ni_s*t(one_mat))%*%exp(eta_next)
    grad_vec <- - eta_next/(sigma_eta^2) - term_exp_eta + 
                      term_y + term_mu_sigma
    eta_hessian <- - (1/(sigma_eta^2))%*%identity_mat - diag(term_exp_eta)
    eta_initial <- eta_next
  }
  mu_grad_vec <- mu_initial*I*(1/(sigma_eta^2) + 1/(c^2)) - 
                           sum(abs(eta_initial))/(sum(sigma_eta^2))
  mu_hessian <- I*(1/(sigma_eta^2) + 1/(c^2))
  while (sum(mu_grad_vec^2) > 0.0001) 
    {
    mu_next <- mu_initial - mu_grad_vec/mu_hessian
    mu_grad_vec <- mu_next*I*(1/(sigma_eta^2) + 1/(c^2)) - 
      sum(abs(eta_initial))/(sum(sigma_eta^2))
    mu_initial <- mu_next
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
