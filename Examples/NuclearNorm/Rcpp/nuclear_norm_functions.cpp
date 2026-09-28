// Compile via nuclear_norm_functions.R, using Rcpp::sourceCpp.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <string>

namespace {
using arma::vec;

arma::uword square_dimension(const vec& x) {
  const arma::uword n = static_cast<arma::uword>(std::sqrt(x.n_elem));
  if (n == 0 || n * n != x.n_elem || !x.is_finite())
    Rcpp::stop("Input must be a finite, nonempty vector of square length.");
  return n;
}

void positive(double value, const char* name) {
  if (!std::isfinite(value) || value <= 0)
    Rcpp::stop("%s must be finite and positive.", name);
}

void validate(const vec& x, const vec& y, double sigma2, double alpha) {
  square_dimension(x);
  if (y.n_elem != x.n_elem || !y.is_finite())
    Rcpp::stop("x/start and y must have equal lengths and finite values.");
  positive(sigma2, "sigma2");
  if (!std::isfinite(alpha) || alpha < 0)
    Rcpp::stop("alpha must be finite and nonnegative.");
}

double nuclear_norm(const vec& x, arma::uword n) {
  vec d;
  if (!arma::svd(d, arma::reshape(x, n, n))) Rcpp::stop("SVD failed.");
  return arma::accu(d);
}

double target(const vec& x, const vec& y, double sigma2, double alpha,
              arma::uword n) {
  return -alpha * nuclear_norm(x, n) - arma::accu(arma::square(y - x)) / (2 * sigma2);
}

struct Prox {
  vec point;
  double norm;
};

Prox proximal(const vec& x, const vec& y, double lambda, double sigma2,
              double alpha, arma::uword n) {
  arma::mat u, v;
  vec d;
  const arma::mat a = arma::reshape((lambda * y + sigma2 * x) / (lambda + sigma2), n, n);
  if (!a.is_finite() || !arma::svd(u, d, v, a))
    Rcpp::stop("Proximity SVD failed; check the step size and inputs.");
  d = arma::clamp(d - alpha * sigma2 * lambda / (lambda + sigma2),
                  0.0, arma::datum::inf);
  return {arma::vectorise(u * arma::diagmat(d) * v.t()), arma::accu(d)};
}

double envelope(const Prox& p, const vec& x, const vec& y,
                double lambda, double sigma2, double alpha) {
  return -alpha * p.norm - arma::accu(arma::square(p.point - x)) / (2 * lambda)
    - arma::accu(arma::square(y - p.point)) / (2 * sigma2);
}

vec normals(arma::uword size) {
  vec out(size);
  for (arma::uword j = 0; j < size; ++j) out[j] = R::rnorm(0, 1);
  return out;
}

double logistic(double x) {
  if (x >= 0) return 1 / (1 + std::exp(-x));
  const double e = std::exp(x);
  return e / (1 + e);
}

double softplus(double x) {
  return std::max(x, 0.0) + std::log1p(std::exp(-std::abs(x)));
}

vec barker_proposal(const vec& x, const vec& gradient, double delta) {
  const vec z = std::sqrt(delta) * normals(x.n_elem);
  vec out(x.n_elem);
  // Draw all normals before the uniforms, matching the original R code.
  for (arma::uword j = 0; j < x.n_elem; ++j)
    out[j] = x[j] + (R::runif(0, 1) <= logistic(z[j] * gradient[j]) ? z[j] : -z[j]);
  return out;
}

double normal_log_density(const vec& x, const vec& mean, double delta) {
  double ans = 0;
  for (arma::uword j = 0; j < x.n_elem; ++j)
    ans += R::dnorm(x[j], mean[j], std::sqrt(delta), true);
  return ans;
}

double barker_log_density(const vec& x, const vec& prop, const vec& gradient,
                          double delta) {
  double ans = 0;
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    const double diff = prop[j] - x[j];
    // As in R, omit the constant length(x)*log(2), which cancels in MH.
    ans += R::dnorm(diff, 0, std::sqrt(delta), true) - softplus(-gradient[j] * diff);
  }
  return ans;
}
} // namespace

// One shared engine preserves the MY/Px target distinction for all proposals.
// [[Rcpp::export]]
SEXP nn_sample_cpp(const arma::vec& y, double alpha, double lambda, double sigma2,
                   int iter, double step, int L, const arma::vec& start,
                   std::string method, bool moreau, bool verbose = true) {
  validate(start, y, sigma2, alpha);
  positive(lambda, "lambda");
  positive(step, "delta/eps_hmc");
  if (iter < 1 || L < 1) Rcpp::stop("iter and L must be positive integers.");
  if (method != "mala" && method != "barker" && method != "hmc")
    Rcpp::stop("Unknown sampler method.");
  const arma::uword n = square_dimension(start), d = start.n_elem;
  Rcpp::NumericMatrix samples(iter, d);
  Rcpp::NumericVector log_weights(moreau ? iter : 0);
  vec current = start;
  Prox prox = proximal(current, y, lambda, sigma2, alpha, n);
  vec gradient = (prox.point - current) / lambda;
  double log_target = moreau ? envelope(prox, current, y, lambda, sigma2, alpha)
    : target(current, y, sigma2, alpha, n);
  double log_weight = moreau ? target(current, y, sigma2, alpha, n) - log_target : 0;
  for (arma::uword j = 0; j < d; ++j) samples(0, j) = current[j];
  if (moreau) log_weights[0] = log_weight;

  const std::string sampler_name = std::string(moreau ? "MY" : "Px")
    + (method == "mala" ? "MALA" : method == "barker" ? "Barker" : "HMC");
  const int progress_updates = std::min(10, iter);
  int next_progress = 1;
  auto report_progress = [&](int completed) {
    // Ten evenly spaced milestones, including 100%, even when iter is not
    // divisible by ten. For iter < 10, report each stored state once.
    if (verbose && next_progress <= progress_updates &&
        static_cast<long long>(completed) * progress_updates >=
        static_cast<long long>(next_progress) * iter) {
      Rcpp::Rcout << sampler_name << ": " << 100 * next_progress / progress_updates
                  << "% (iteration " << completed << " / " << iter << ")\n";
      Rcpp::Rcout.flush();
      ++next_progress;
    }
  };
  report_progress(1);

  // Preserve R's column-major, upfront HMC momentum draws for seeded parity.
  arma::mat momenta;
  if (method == "hmc") {
    momenta.set_size(iter, d);
    for (arma::uword j = 0; j < d; ++j) {
      Rcpp::checkUserInterrupt();
      for (int i = 0; i < iter; ++i) momenta(i, j) = R::rnorm(0, 1);
    }
  }
  int accepted = 0;
  for (int i = 1; i < iter; ++i) {
    Rcpp::checkUserInterrupt();
    vec proposal, next_gradient;
    Prox next_prox;
    double correction = 0;
    if (method == "hmc") {
      const vec initial_momentum = momenta.row(i).t();
      vec momentum = initial_momentum + step * gradient / 2;
      proposal = current;
      const int steps = R::runif(0, 1) <= 0.05 ? 1 : L;
      for (int j = 0; j < steps; ++j) {
        Rcpp::checkUserInterrupt();
        proposal += step * momentum;
        next_prox = proximal(proposal, y, lambda, sigma2, alpha, n);
        next_gradient = (next_prox.point - proposal) / lambda;
        if (j + 1 != steps) momentum += step * next_gradient;
      }
      momentum += step * next_gradient / 2;
      momentum = -momentum;
      correction = (arma::dot(initial_momentum, initial_momentum)
                    - arma::dot(momentum, momentum)) / 2;
    } else {
      vec mean;
      if (method == "mala") {
        mean = current + step * gradient / 2;
        proposal = mean + std::sqrt(step) * normals(d);
      } else {
        proposal = barker_proposal(current, gradient, step);
      }
      next_prox = proximal(proposal, y, lambda, sigma2, alpha, n);
      next_gradient = (next_prox.point - proposal) / lambda;
      correction = method == "mala"
        ? normal_log_density(current, proposal + step * next_gradient / 2, step)
          - normal_log_density(proposal, mean, step)
        : barker_log_density(proposal, current, next_gradient, step)
          - barker_log_density(current, proposal, gradient, step);
    }
    const double next_target = moreau
      ? envelope(next_prox, proposal, y, lambda, sigma2, alpha)
      : target(proposal, y, sigma2, alpha, n);
    if (std::log(R::runif(0, 1)) <= next_target - log_target + correction) {
      current = std::move(proposal);
      prox = std::move(next_prox);
      gradient = std::move(next_gradient);
      log_target = next_target;
      if (moreau) log_weight = target(current, y, sigma2, alpha, n) - log_target;
      ++accepted;
    }
    for (arma::uword j = 0; j < d; ++j) samples(i, j) = current[j];
    if (moreau) log_weights[i] = log_weight;
    report_progress(i + 1);
  }
  // Keep the original accept/iter convention and positional return formats.
  const double rate = static_cast<double>(accepted) / iter;
  if (verbose) {
    Rcpp::Rcout << sampler_name << " acceptance rate: ";
    if (iter > 1)
      Rcpp::Rcout << static_cast<double>(accepted) / (iter - 1)
                  << " (" << accepted << " / " << (iter - 1) << " proposals accepted)\n";
    else
      Rcpp::Rcout << "N/A (no proposals; only the starting point was stored)\n";
    Rcpp::Rcout.flush();
  }
  if (moreau && method == "mala") return Rcpp::List::create(samples, log_weights);
  if (moreau) return Rcpp::List::create(samples, log_weights, rate);
  if (method == "mala") return samples;
  return Rcpp::List::create(samples, rate);
}
