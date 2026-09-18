# -------------------------------------------------------------------------
# Cox-based semiparametric quantile regression estimators
# -------------------------------------------------------------------------
# Includes:
#   1) No covariates
#   2) One additional covariate (X2)
#
# Required package: survival
# -------------------------------------------------------------------------

library(survival)

# -------------------------------------------------------------------------
# 1) No covariates
# -------------------------------------------------------------------------
# Fits: Surv(event_time, status) ~ treatment
# Returns:
#   - beta_hat(tau) = Q_tau(Y | treatment = 1) - Q_tau(Y | treatment = 0)
#   - estimated quantiles for treatment = 0 and 1

cox_qr_estimate <- function(data, tau) {
  
  required_cols <- c("event_time", "status", "treatment")
  if (!all(required_cols %in% names(data))) {
    stop("Data must contain columns: event_time, status, treatment.")
  }
  
  tau <- as.numeric(tau)
  if (any(is.na(tau)) || any(tau <= 0 | tau >= 1)) {
    stop("Values in 'tau' must be strictly between 0 and 1.")
  }
  
  if (!all(stats::na.omit(data$treatment) %in% c(0, 1))) {
    stop("'treatment' must be coded as 0/1.")
  }
  
  cox_model <- coxph(Surv(event_time, status) ~ treatment, data = data, ties = "efron")
  
  base_hazard <- basehaz(cox_model, centered = FALSE)
  base_hazard <- base_hazard[order(base_hazard$time), , drop = FALSE]
  
  beta_treatment <- unname(coef(cox_model)["treatment"])
  
  invert_baseline_hazard <- function(u) {
    approx(
      x = base_hazard$hazard,
      y = base_hazard$time,
      xout = u,
      method = "linear",
      rule = 2
    )$y
  }
  
  H0 <- log(1 / (1 - tau))
  
  q0 <- invert_baseline_hazard(H0)
  q1 <- invert_baseline_hazard(H0 / exp(beta_treatment))
  
  beta_hat <- q1 - q0
  
  list(
    coef_values = data.frame(
      tau = tau,
      beta_hat = beta_hat
    ),
    inv_values = data.frame(
      tau = tau,
      inv_treatment_0 = q0,
      inv_treatment_1 = q1
    ),
    cox_model = cox_model,
    base_hazard = base_hazard
  )
}

# -------------------------------------------------------------------------
# 2) One additional covariate
# -------------------------------------------------------------------------
# Fits: Surv(event_time, status) ~ treatment + X2
# Computes:
#   beta_hat(tau, x2_i) = Q_tau(Y | treatment = 1, X2 = x2_i) -
#                         Q_tau(Y | treatment = 0, X2 = x2_i)
# Returns:
#   - beta_hat(tau, x2_i) for each individual
#   - AQE_hat(tau) = sample average of beta_hat(tau, x2_i)

cox_qr_estimate_cov <- function(data, tau) {
  
  required_cols <- c("event_time", "status", "treatment", "X2")
  if (!all(required_cols %in% names(data))) {
    stop("Data must contain columns: event_time, status, treatment, X2.")
  }
  
  tau <- as.numeric(tau)
  if (any(is.na(tau)) || any(tau <= 0 | tau >= 1)) {
    stop("Values in 'tau' must be strictly between 0 and 1.")
  }
  
  if (!all(stats::na.omit(data$treatment) %in% c(0, 1))) {
    stop("'treatment' must be coded as 0/1.")
  }
  
  cox_model <- coxph(Surv(event_time, status) ~ treatment + X2, data = data, ties = "efron")
  
  base_hazard <- basehaz(cox_model, centered = FALSE)
  base_hazard <- base_hazard[order(base_hazard$time), , drop = FALSE]
  
  beta <- coef(cox_model)
  beta_treatment <- unname(beta["treatment"])
  beta_X2 <- unname(beta["X2"])
  
  invert_baseline_hazard <- function(u) {
    approx(
      x = base_hazard$hazard,
      y = base_hazard$time,
      xout = u,
      method = "linear",
      rule = 2
    )$y
  }
  
  n <- nrow(data)
  
  out <- expand.grid(
    tau = tau,
    id = seq_len(n)
  )
  out <- out[order(out$tau, out$id), , drop = FALSE]
  out$X2 <- data$X2[out$id]
  
  log_term <- log(1 / (1 - out$tau))
  eta_X2 <- beta_X2 * out$X2
  
  haz0 <- log_term / exp(eta_X2)
  haz1 <- log_term / exp(beta_treatment + eta_X2)
  
  hmin <- min(base_hazard$hazard, na.rm = TRUE)
  hmax <- max(base_hazard$hazard, na.rm = TRUE)
  extrapolation_count <- sum(
    haz0 < hmin | haz0 > hmax | haz1 < hmin | haz1 > hmax
  )
  
  q0 <- invert_baseline_hazard(haz0)
  q1 <- invert_baseline_hazard(haz1)
  
  out$inv_treatment_0 <- q0
  out$inv_treatment_1 <- q1
  out$beta_hat <- q1 - q0
  
  aqe_hat <- aggregate(beta_hat ~ tau, data = out, FUN = mean)
  names(aqe_hat)[2] <- "AQE_hat"
  
  list(
    coef_values = out,
    aqe_hat = aqe_hat,
    cox_model = cox_model,
    base_hazard = base_hazard,
    extrapolation_count = extrapolation_count
  )
}
