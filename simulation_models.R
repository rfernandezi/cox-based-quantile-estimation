# -------------------------------------------------------------------------
# Data-generating models for Monte Carlo simulations
# -------------------------------------------------------------------------
# Includes:
#   1) Models I-VI: no covariates
#   2) Models VII-VIII: one additional covariate (X2)
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Models I-VI (no covariates)
# -------------------------------------------------------------------------
# Each model defines two quantile functions:
#   - f0(tau): quantile function under treatment = 0
#   - f1(tau): quantile function under treatment = 1

get_models_I_VI <- function() {
  list(
    list(
      name = "Model I",
      f0 = function(t) qexp(t, rate = 1),
      f1 = function(t) qexp(t, rate = 1)
    ),
    list(
      name = "Model II",
      f0 = function(t) qexp(t, rate = 1),
      f1 = function(t) qexp(t, rate = 1 / 2)
    ),
    list(
      name = "Model III",
      f0 = function(t) qexp(t, rate = 1),
      f1 = function(t) qnorm(t, mean = 2, sd = 3 / 2)
    ),
    list(
      name = "Model IV",
      f0 = function(t) qnorm(t, mean = 8, sd = 1),
      f1 = function(t) qnorm(t, mean = 8, sd = 2)
    ),
    list(
      name = "Model V",
      f0 = function(t) qnorm(t, mean = 0, sd = 1),
      f1 = function(t) qnorm(t, mean = 1, sd = 1)
    ),
    list(
      name = "Model VI",
      f0 = function(t) qweibull(t, shape = 1 / 2, scale = 1 / 2),
      f1 = function(t) qweibull(t, shape = 3 / 4, scale = 1)
    )
  )
}


# -------------------------------------------------------------------------
# Helper functions for Models VII-VIII
# -------------------------------------------------------------------------

generate_cox_times_cov <- function(u, lambda, beta1, beta2, X1, X2) {
  -log(1 - u) / (lambda * exp(beta1 * X1 + beta2 * X2))
}

generate_qr_times <- function(u, X1, X2) {
  qnorm(u, mean = 2, sd = 1) +
    qnorm(u, mean = 0.5, sd = 0.2) * X1 +
    qnorm(u, mean = 1, sd = 0.3) * X2
}

get_aqe_closed <- function(tau, mu0, sigma0, mu1, sigma1, pi, beta1, beta2, lambda) {
  constant <- (-log(1 - tau) / lambda) * (exp(-beta1) - 1)
  mixture_term <- (1 - pi) * exp(-beta2 * mu0 + 0.5 * beta2^2 * sigma0^2) +
    pi * exp(-beta2 * mu1 + 0.5 * beta2^2 * sigma1^2)
  constant * mixture_term
}


# -------------------------------------------------------------------------
# Models VII-VIII (one additional covariate)
# -------------------------------------------------------------------------
# Returns:
#   - data-generating mechanism
#   - true AQE(tau)

get_models_VII_VIII <- function(taus, p, lambda = 2, beta1 = 1, beta2 = 0.5) {
  
  list(
    list(
      name = "Model VII-A",
      generate_X2 = function(X1) rnorm(length(X1), mean = 0, sd = 1),
      generate_times = function(u, X1, X2) {
        generate_cox_times_cov(u, lambda, beta1, beta2, X1, X2)
      },
      true_aqe = sapply(
        taus,
        function(t) get_aqe_closed(
          tau = t,
          mu0 = 0, sigma0 = 1,
          mu1 = 0, sigma1 = 1,
          pi = p, beta1 = beta1, beta2 = beta2, lambda = lambda
        )
      )
    ),
    
    list(
      name = "Model VII-B",
      generate_X2 = function(X1) ifelse(
        X1 == 0,
        rnorm(length(X1), mean = 0, sd = 1),
        rnorm(length(X1), mean = 1, sd = 1)
      ),
      generate_times = function(u, X1, X2) {
        generate_cox_times_cov(u, lambda, beta1, beta2, X1, X2)
      },
      true_aqe = sapply(
        taus,
        function(t) get_aqe_closed(
          tau = t,
          mu0 = 0, sigma0 = 1,
          mu1 = 1, sigma1 = 1,
          pi = p, beta1 = beta1, beta2 = beta2, lambda = lambda
        )
      )
    ),
    
    list(
      name = "Model VII-C",
      generate_X2 = function(X1) ifelse(
        X1 == 0,
        rnorm(length(X1), mean = 0, sd = 1),
        rnorm(length(X1), mean = 0, sd = 3 / 2)
      ),
      generate_times = function(u, X1, X2) {
        generate_cox_times_cov(u, lambda, beta1, beta2, X1, X2)
      },
      true_aqe = sapply(
        taus,
        function(t) get_aqe_closed(
          tau = t,
          mu0 = 0, sigma0 = 1,
          mu1 = 0, sigma1 = 3 / 2,
          pi = p, beta1 = beta1, beta2 = beta2, lambda = lambda
        )
      )
    ),
    
    list(
      name = "Model VIII-A",
      generate_X2 = function(X1) rnorm(length(X1), mean = 0, sd = 1),
      generate_times = function(u, X1, X2) {
        generate_qr_times(u, X1, X2)
      },
      true_aqe = qnorm(taus, mean = 0.5, sd = 0.2)
    ),
    
    list(
      name = "Model VIII-B",
      generate_X2 = function(X1) ifelse(
        X1 == 0,
        rnorm(length(X1), mean = 0, sd = 1),
        rnorm(length(X1), mean = 1, sd = 1)
      ),
      generate_times = function(u, X1, X2) {
        generate_qr_times(u, X1, X2)
      },
      true_aqe = qnorm(taus, mean = 0.5, sd = 0.2)
    ),
    
    list(
      name = "Model VIII-C",
      generate_X2 = function(X1) ifelse(
        X1 == 0,
        rnorm(length(X1), mean = 0, sd = 1),
        rnorm(length(X1), mean = 0, sd = 3 / 2)
      ),
      generate_times = function(u, X1, X2) {
        generate_qr_times(u, X1, X2)
      },
      true_aqe = qnorm(taus, mean = 0.5, sd = 0.2)
    )
  )
}