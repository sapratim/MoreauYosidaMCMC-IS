####  Poisson random effects model
#### The function evaluated below of targets and the derivatives are omitting constant
#### of proportionality.

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
tol_nr <- 1e-5
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
  print(iter)
  optima <- c(eta_next, mu_next)
  return(optima)
}

MY_env <- function(eta, mu, data, optima, lambda)
{
  term <- rbind(eta, mu)
  value <- true_target(eta, mu, data) + sum(term^2)/lambda
  return(-value)
}

eta_start <- eta_vec + rnorm(I, 0, 2)
mu_start <- 5
check <- proxfunc(eta_vec, mu, lambda = lambda, eta_initial = eta_start, 
         mu_initial = mu_start, sigma = sigma_eta)

cbind(check, c(eta_vec,mu))
