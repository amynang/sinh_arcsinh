# This is the Sinh-Arcsinh distribution (Jones & Pewsey, 2009; 10.1093/biomet/asp053)
# based on Tim Wolock's definition as a custom family for brms
# Wolock et al. 2021
# https://elifesciences.org/articles/68318
# The original can be found here:
# https://github.com/twolock/distreg-illustration/blob/main/R/stan_funs.R

# I only made minor changes in syntax following 
# https://github.com/paul-buerkner/custom-brms-families/tree/master

# Defining a density function and rng in R allows us to bypass expose_functions()
# which can create trouble in Linux (and crashes RStudio in my case)

# density
dsinhasinh <- function(y, mu, sigma, epsilon, delta, log = FALSE) {
  nll = 0
  sigma_star = sigma * delta
  y_z = (y - mu)/sigma_star
  
  S_y = sinh(epsilon + delta * asinh(y_z))
  S_y_2 = S_y * S_y
  C_y = sqrt(1 + S_y_2)
  nll = nll -0.5 * S_y_2 - log(sigma_star)
  nll = nll + log(delta) + log(C_y) - log(sqrt(1 + y_z*y_z))
  
  if (log) {
    return(nll)
  } else {
    return(exp(nll))
  }
}

# rng
rsinhasinh <- function(n, mu, sigma, epsilon, delta) {
  mu + sigma * delta * sinh((asinh(rnorm(n, 0, 1)) - epsilon)/delta)
}

log_lik_sinhasinh <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  epsilon <- brms::get_dpar(prep, "epsilon", i = i)
  delta <- brms::get_dpar(prep, "delta", i = i)
  y <- prep$data$Y[i]
  #sinhasinh_lpdf(y, mu, sigma, epsilon, delta)
  return(dsinhasinh(y = y, mu = mu, sigma = sigma, epsilon = epsilon, delta = delta, log = TRUE))
}

posterior_predict_sinhasinh <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  epsilon <- brms::get_dpar(prep, "epsilon", i = i)
  delta <- brms::get_dpar(prep, "delta", i = i)
  #sinhasinh_rng(mu, sigma, epsilon, delta)
  return(rsinhasinh(n = prep$ndraws, mu = mu, sigma = sigma, epsilon = epsilon, delta = delta))
}

posterior_epred_sinhasinh <- function(prep) {
  brms::get_dpar(prep, "mu")
}

sinhasinh <- function(link = "identity", 
                      link_sigma = "log",
                      link_epsilon = "identity",
                      link_delta = "log") {
  custom_family(
    "sinhasinh",
    dpars = c("mu", "sigma", "epsilon", "delta"),
    links = c(link, link_sigma, link_epsilon, link_delta),
    lb = c(NA, NA, NA, NA),
    ub = c(NA, NA, NA, NA),
    type = "real",
    log_lik = log_lik_sinhasinh,
    posterior_predict = posterior_predict_sinhasinh,
    posterior_epred = posterior_epred_sinhasinh,
  )
}

stan_funs <- "
  real sinhasinh_lpdf(real y, real mu, real sigma, real epsilon, real delta) {
    real y_z;
    real sigma_star;
    real S_y;
    real S_y_2;
    real C_y;
    real nll;
    
    nll = 0;
    sigma_star = sigma * delta;
    y_z = (y - mu)/sigma_star;
    
    S_y = sinh(epsilon + delta * asinh(y_z));
    S_y_2 = S_y * S_y;
    C_y = sqrt(1 + S_y_2);
    nll += -0.5 * S_y_2 - log(sigma_star);
    nll += log(delta) + log(C_y) - log(sqrt(1 + y_z*y_z));
    return nll;
  }
  real sinhasinh_rng(real mu, real sigma, real epsilon, real delta) {
    return (mu + sigma * delta * sinh((asinh(normal_rng(0, 1)) - epsilon)/delta));
  }
"