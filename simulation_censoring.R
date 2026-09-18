# =============================================================================
# Simulation study with right censoring
# =============================================================================
#
# Sensitivity analysis evaluating the performance of Cox-QR under
# 0%, 20%, and 40% right censoring.
#
# Models II and VII-C satisfy the proportional hazards assumption,
# whereas Models VI and VIII-C violate it.
#
# Cox-QR is compared with censored quantile regression using the
# Peng-Huang and Portnoy methods implemented in the quantreg package.
# =============================================================================

library(survival)
library(quantreg)

source("cox_qr_estimator.R")
source("simulation_models.R")

set.seed(123)

taus <- c(0.25, 0.50, 0.75)


# -----------------------------------------------------------------------------
# Right-censoring mechanism
# -----------------------------------------------------------------------------

# Calibrate the rate of an exponential censoring distribution to obtain
# approximately the desired censoring proportion.
calc_lambda_c <- function(Tcal, target_censoring) {
  
  f <- function(lambda_c) {
    mean(
      (1 - exp(-lambda_c * Tcal)) * (Tcal > 0)
    ) - target_censoring
  }
  
  uniroot(
    f,
    interval = c(1e-10, 100)
  )$root
}


# Apply independent right censoring:
# C ~ Exponential(lambda_c)
# observed time = min(T, C)
# status = 1(T <= C)
apply_censoring <- function(event_time, lambda_c) {
  
  censor_time <- rexp(
    length(event_time),
    rate = lambda_c
  )
  
  observed_time <- pmin(
    event_time,
    censor_time
  )
  
  status <- as.integer(
    event_time <= censor_time
  )
  
  data.frame(
    event_time = observed_time,
    status = status
  )
}


# -----------------------------------------------------------------------------
# Censored quantile regression comparators
# -----------------------------------------------------------------------------

fit_censored_qr <- function(data, taus, method = "PengHuang",
                            covariate = FALSE) {
  
  if (covariate) {
    
    fit <- crq(
      Surv(event_time, status) ~ treatment + X2,
      data = data,
      method = method
    )
    
  } else {
    
    fit <- crq(
      Surv(event_time, status) ~ treatment,
      data = data,
      method = method
    )
  }
  
  sapply(taus, function(tau) {
    
    cf <- tryCatch(
      coef(fit, taus = tau),
      error = function(e) NULL
    )
    
    if (is.null(cf)) {
      return(NA_real_)
    }
    
    if (is.matrix(cf)) {
      cf <- cf[, 1]
    }
    
    if (!"treatment" %in% names(cf)) {
      return(NA_real_)
    }
    
    as.numeric(cf["treatment"])
  })
}


# -----------------------------------------------------------------------------
# Helper: generate censored data
# -----------------------------------------------------------------------------

add_censoring <- function(data, target_censoring) {
  
  # No censoring
  if (target_censoring == 0) {
    data$status <- 1
    return(data)
  }
  
  lambda_c <- calc_lambda_c(
    data$event_time,
    target_censoring
  )
  
  cens <- apply_censoring(
    data$event_time,
    lambda_c
  )
  
  data$event_time <- cens$event_time
  data$status <- cens$status
  
  data
}


# -----------------------------------------------------------------------------
# Simulation settings
# -----------------------------------------------------------------------------

censoring_levels <- c(0, 0.20, 0.40)

sample_configs <- list(
  size1 = list(n = 100, p = 1 / 2),
  size2 = list(n = 300, p = 1 / 2),
  size3 = list(n = 300, p = 1 / 3)
)

n_sim <- 2000
