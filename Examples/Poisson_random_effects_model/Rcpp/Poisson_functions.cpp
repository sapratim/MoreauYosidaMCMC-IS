// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {
using Vec = std::vector<double>;
void positive(double x, const char* name) {
  if (!std::isfinite(x) || x <= 0) Rcpp::stop("%s must be finite and positive.", name);
}
struct Model {
  int n, m;
  double a, b;
  Vec sums;
  Model(Rcpp::NumericMatrix data, double sigma, double prior_sd) :
    n(data.nrow()), m(data.ncol()), sums(n, 0) {
    positive(sigma, "sigma"); positive(prior_sd, "prior_sd");
    a = 1 / (sigma * sigma); b = 1 / (prior_sd * prior_sd);
    positive(a, "inverse sigma squared"); positive(b, "inverse prior variance");
    if (!n || !m) Rcpp::stop("data must be a nonempty matrix.");
    for (int j = 0; j < m; ++j) for (int i = 0; i < n; ++i) {
      double y = data(i, j);
      if (!std::isfinite(y) || y < 0 || y != std::floor(y))
        Rcpp::stop("data must contain finite nonnegative counts.");
      sums[i] += y;
      if (!std::isfinite(sums[i])) Rcpp::stop("Data row sums overflow.");
    }
  }
  void validate(const Vec& x) const {
    if (x.size() != static_cast<size_t>(n + 1)) Rcpp::stop("State length must be nrow(data) + 1.");
    for (double z : x) if (!std::isfinite(z)) Rcpp::stop("State must be finite.");
  }
  double target(const Vec& x) const {
    double value = -0.5 * b * x[n] * x[n];
    for (int j = 0; j < n; ++j) {
      double diff = x[j] - x[n];
      value += sums[j] * x[j] - m * std::exp(x[j]) - 0.5 * a * diff * diff;
    }
    return value;
  }
  Vec gradient(const Vec& x) const {
    Vec g(n + 1); g[n] = -b * x[n];
    for (int j = 0; j < n; ++j) {
      double diff = a * (x[j] - x[n]);
      g[j] = sums[j] - m * std::exp(x[j]) - diff;
      g[n] += diff;
    }
    return g;
  }
};
struct Prox { Vec point, rates; };
// Full Newton solve of the arrowhead Hessian: O(n) work, no dense matrices.
// Initialization depends only on x, so the numerical force is deterministic.
Prox proximal(const Model& model, const Vec& x, double lambda, double tol) {
  const int n = model.n;
  const double il = 1 / lambda, a = model.a;
  Prox out{ x, Vec(n) };
  Vec& p = out.point;
  double sum = 0;
  for (int j = 0; j < n; ++j) {
    p[j] = std::min(x[j], std::log((model.sums[j] + 1) / model.m));
    sum += p[j];
  }
  const double hmu = n * a + model.b + il;
  p[n] = (a * sum + il * x[n]) / hmu;
  Vec g(n + 1), diag(n), step(n + 1);
  for (int k = 0; k < 100; ++k) {
    Rcpp::checkUserInterrupt();
    double sg = 0, si = 0, norm = 0, scale = 1;
    g[n] = model.b * p[n] + il * (p[n] - x[n]);
    for (int j = 0; j < n; ++j) {
      out.rates[j] = model.m * std::exp(p[j]);
      g[j] = a * (p[j] - p[n]) + out.rates[j] - model.sums[j] + il * (p[j] - x[j]);
      g[n] += a * (p[n] - p[j]);
      diag[j] = a + out.rates[j] + il;
      sg += g[j] / diag[j]; si += 1 / diag[j];
      norm = std::max(norm, std::abs(g[j]));
      scale = std::max(scale, out.rates[j] + model.sums[j] +
        a * (std::abs(p[j]) + std::abs(p[n])) + il * (std::abs(p[j]) + std::abs(x[j])));
    }
    norm = std::max(norm, std::abs(g[n]));
    scale = std::max(scale, hmu * std::abs(p[n]) + a * std::abs(sum) + il * std::abs(x[n]));
    if (norm <= tol + 32 * std::numeric_limits<double>::epsilon() * scale) return out;
    step[n] = -(g[n] + a * sg) / (hmu - a * a * si);
    double slope = g[n] * step[n];
    for (int j = 0; j < n; ++j) {
      step[j] = (-g[j] + a * step[n]) / diag[j];
      slope += g[j] * step[j];
    }
    if (!std::isfinite(slope) || slope >= 0) Rcpp::stop("Proximal Newton solve failed.");
    double t = 1;
    bool accepted = false;
    for (int back = 0; back < 60; ++back, t *= 0.5) {
      // Compute objective changes directly to avoid subtracting large targets.
      const double dm = t * step[n];
      double change = model.b * dm * (p[n] + 0.5 * dm) +
        il * dm * (p[n] - x[n] + 0.5 * dm);
      for (int j = 0; j < n; ++j) {
        const double de = t * step[j], dr = de - dm;
        change += out.rates[j] * std::expm1(de) - model.sums[j] * de +
          a * dr * (p[j] - p[n] + 0.5 * dr) + il * de * (p[j] - x[j] + 0.5 * de);
      }
      if (std::isfinite(change) && change <= 1e-4 * t * slope) { accepted = true; break; }
    }
    if (!accepted) Rcpp::stop("Proximal line search failed to converge.");
    sum = 0;
    for (int j = 0; j <= n; ++j) { p[j] += t * step[j]; if (j < n) sum += p[j]; }
  }
  Rcpp::stop("Proximal solver exceeded 100 Newton iterations.");
  return out;
}
double envelope(const Model& m, const Vec& x, const Prox& p, double lambda) {
  double value = -0.5 * m.b * p.point[m.n] * p.point[m.n];
  for (int j = 0; j <= m.n; ++j) {
    const double d = p.point[j] - x[j]; value -= d * d / (2 * lambda);
    if (j < m.n) {
      const double r = p.point[j] - p.point[m.n];
      value += m.sums[j] * p.point[j] - p.rates[j] - 0.5 * m.a * r * r;
    }
  }
  return value;
}
// log p(x) - log p_lambda(x), using cached prox exponentials and expm1.
// Only called for the initial state and accepted MY proposals.
double weight(const Model& m, const Vec& x, const Prox& p, double lambda) {
  const double dm = x[m.n] - p.point[m.n];
  double value = -m.b * dm * (p.point[m.n] + 0.5 * dm) + dm * dm / (2 * lambda);
  for (int j = 0; j < m.n; ++j) {
    const double de = x[j] - p.point[j], dr = de - dm;
    value += m.sums[j] * de - p.rates[j] * std::expm1(de) -
      m.a * dr * (p.point[j] - p.point[m.n] + 0.5 * dr) + de * de / (2 * lambda);
  }
  return value;
}
Vec force(const Vec& x, const Prox& p, double lambda) {
  Vec g(x.size());
  for (size_t j = 0; j < x.size(); ++j) g[j] = (p.point[j] - x[j]) / lambda;
  return g;
}
double softplus(double x) { return std::max(x, 0.0) + std::log1p(std::exp(-std::abs(x))); }
double logistic(double x) {
  if (x >= 0) return 1 / (1 + std::exp(-x));
  double e = std::exp(x); return e / (1 + e);
}
}

// [[Rcpp::export]]
Rcpp::NumericVector poisson_prox_cpp(Rcpp::NumericVector state, Rcpp::NumericMatrix data,
    double lambda, double sigma, double prior_sd, double tol = 1e-8) {
  positive(lambda, "lambda"); positive(1 / lambda, "inverse lambda"); positive(tol, "tol");
  Model m(data, sigma, prior_sd); Vec x = Rcpp::as<Vec>(state); m.validate(x);
  return Rcpp::wrap(proximal(m, x, lambda, tol).point);
}

// All seven public R samplers dispatch to this compiled engine.
// [[Rcpp::export]]
SEXP poisson_sample_cpp(Rcpp::NumericVector start, Rcpp::NumericMatrix data,
    double lambda, double sigma, double prior_sd, int iter, double step, int L,
    std::string method, bool moreau, bool true_gradient = false,
    bool verbose = true, double tol = 1e-8) {
  Model m(data, sigma, prior_sd); Vec current = Rcpp::as<Vec>(start); m.validate(current);
  positive(lambda, "lambda"); positive(1 / lambda, "inverse lambda"); positive(step, "step"); positive(tol, "tol");
  if (iter < 1 || L < 1) Rcpp::stop("iter and L must be positive integers.");
  if (method != "mala" && method != "barker" && method != "hmc") Rcpp::stop("Unknown method.");
  if (true_gradient && (moreau || method != "barker")) Rcpp::stop("True gradient is supported for exact Barker only.");
  const int d = m.n + 1;
  Rcpp::NumericMatrix samples(iter, d);
  Rcpp::NumericVector weights(moreau ? iter : 0);
  Prox prox;
  if (!true_gradient) prox = proximal(m, current, lambda, tol);
  Vec gradient = true_gradient ? m.gradient(current) : force(current, prox, lambda);
  double target = moreau ? envelope(m, current, prox, lambda) : m.target(current);
  if (!std::isfinite(target)) Rcpp::stop("Initial log target is not finite.");
  double log_weight = moreau ? weight(m, current, prox, lambda) : 0;
  for (int j = 0; j < d; ++j) samples(0, j) = current[j];
  if (moreau) weights[0] = log_weight;

  const std::string name = std::string(true_gradient ? "" : moreau ? "MY" : "Px") +
    (method == "mala" ? "MALA" : method == "hmc" ? "HMC" : "Barker");
  int milestone = 1, accepted = 0;
  const int updates = std::min(10, iter);
  auto report = [&](int completed) {
    if (verbose && milestone <= updates && static_cast<long long>(completed) * updates >=
        static_cast<long long>(milestone) * iter) {
      Rcpp::Rcout << name << ": " << 100 * milestone / updates << "% (" << completed << " / " << iter << ")\n";
      Rcpp::Rcout.flush(); ++milestone;
    }
  };
  report(1);
  // Match the original R HMC draw order, including the unused first row.
  Rcpp::NumericMatrix momenta(method == "hmc" ? iter : 0, method == "hmc" ? d : 0);
  if (method == "hmc") for (int j = 0; j < d; ++j) {
    Rcpp::checkUserInterrupt();
    for (int i = 0; i < iter; ++i) momenta(i, j) = R::rnorm(0, 1);
  }
  for (int i = 1; i < iter; ++i) {
    Rcpp::checkUserInterrupt();
    Vec proposal(d), next_gradient(d);
    Prox next_prox;
    double correction = 0;
    bool valid = true;
    auto evaluate_force = [&]() {
      for (double v : proposal) if (!std::isfinite(v)) return false;
      if (true_gradient) next_gradient = m.gradient(proposal);
      else { next_prox = proximal(m, proposal, lambda, tol); next_gradient = force(proposal, next_prox, lambda); }
      for (double v : next_gradient) if (!std::isfinite(v)) return false;
      return true;
    };
    if (method == "hmc") {
      Vec momentum(d); proposal = current;
      for (int j = 0; j < d; ++j) {
        correction += momenta(i, j) * momenta(i, j) / 2;
        momentum[j] = momenta(i, j) + step * gradient[j] / 2;
      }
      const int steps = R::runif(0, 1) <= 0.05 ? 1 : L;
      for (int k = 0; k < steps; ++k) {
        Rcpp::checkUserInterrupt();
        for (int j = 0; j < d; ++j) proposal[j] += step * momentum[j];
        if (!(valid = evaluate_force())) break;
        if (k + 1 != steps) for (int j = 0; j < d; ++j) momentum[j] += step * next_gradient[j];
      }
      if (valid) for (int j = 0; j < d; ++j) {
        momentum[j] += step * next_gradient[j] / 2;
        correction -= momentum[j] * momentum[j] / 2;
      }
    } else {
      Vec z(d);
      for (int j = 0; j < d; ++j) z[j] = std::sqrt(step) * R::rnorm(0, 1);
      for (int j = 0; j < d; ++j) proposal[j] = current[j] +
        (method == "mala" ? step * gradient[j] / 2 + z[j] :
          (R::runif(0, 1) <= logistic(z[j] * gradient[j]) ? z[j] : -z[j]));
      valid = evaluate_force();
      if (valid) for (int j = 0; j < d; ++j) {
        const double diff = proposal[j] - current[j];
        if (method == "mala") {
          double reverse = -diff - step * next_gradient[j] / 2;
          double forward = diff - step * gradient[j] / 2;
          correction += (forward * forward - reverse * reverse) / (2 * step);
        } else {
          // Symmetric Gaussian terms and the d*log(2) constants cancel.
          correction += softplus(-gradient[j] * diff) - softplus(next_gradient[j] * diff);
        }
      }
    }
    double next_target = valid ? (moreau ? envelope(m, proposal, next_prox, lambda) : m.target(proposal)) : R_NegInf;
    const double log_u = std::log(R::runif(0, 1));
    if (valid && std::isfinite(next_target) && log_u <= next_target - target + correction) {
      current.swap(proposal); gradient.swap(next_gradient); target = next_target;
      if (moreau) log_weight = weight(m, current, next_prox, lambda);
      ++accepted;
    }
    for (int j = 0; j < d; ++j) samples(i, j) = current[j];
    if (moreau) weights[i] = log_weight;
    report(i + 1);
  }
  if (verbose) {
    Rcpp::Rcout << name << " acceptance rate: ";
    if (iter > 1) Rcpp::Rcout << static_cast<double>(accepted) / (iter - 1) << "\n";
    else Rcpp::Rcout << "N/A (no proposals)\n";
    Rcpp::Rcout.flush();
  }
  // Preserve the original positional outputs and returned accept/iter convention.
  const double rate = static_cast<double>(accepted) / iter;
  if (moreau && method == "mala") return Rcpp::List::create(samples, weights);
  if (moreau) return Rcpp::List::create(samples, weights, rate);
  if (method == "mala") return samples;
  return Rcpp::List::create(samples, rate);
}
