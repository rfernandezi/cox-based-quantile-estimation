# -------------------------------------------------------------------------
# Monte Carlo simulations: coverage of bootstrap confidence intervals
# -------------------------------------------------------------------------
# Computes empirical coverage for:
#   - beta(tau) in Models I-VI
#   - AQE(tau) in Models VII-VIII
# -------------------------------------------------------------------------

library(survival)

source("R/cox_qr_estimator.R")
source("R/simulation_models.R")

set.seed(123)

taus <- c(0.25, 0.50, 0.75)
alpha <- 0.05
B <- 2000


# -------------------------------------------------------------------------
# Helper: bootstrap estimator
# -------------------------------------------------------------------------

bootstrap_cox_qr <- function(data, taus, B) {
  
  n <- nrow(data)
  boot <- matrix(NA, B, length(taus))
  
  for (b in seq_len(B)) {
    
    data_boot <- data[sample(n, replace = TRUE), ]
    
    fit <- cox_qr_estimate(data_boot, taus)
    boot[b, ] <- fit$coef_values$beta_hat
  }
  
  boot
}


bootstrap_cox_qr_cov <- function(data, taus, B) {
  
  n <- nrow(data)
  boot <- matrix(NA, B, length(taus))
  
  for (b in seq_len(B)) {
    
    data_boot <- data[sample(n, replace = TRUE), ]
    
    fit <- cox_qr_estimate_cov(data_boot, taus)
    boot[b, ] <- fit$aqe_hat$AQE_hat
  }
  
  boot
}


# -------------------------------------------------------------------------
# Coverage: Models I-VI
# -------------------------------------------------------------------------

run_coverage_I_VI <- function(n_sim = 2000,
                              n = 300,
                              p = 1/3,
                              taus = taus,
                              out_file = "simulations/coverage_models_I_VI.csv") {
  
  models <- get_models_I_VI()
  results <- vector("list", length(models) * n_sim)
  idx <- 1
  
  for (m in models) {
    
    beta_true <- sapply(taus, function(t) m$f1(t) - m$f0(t))
    
    for (sim in seq_len(n_sim)) {
      
      treatment <- rbinom(n, 1, p)
      event_time <- numeric(n)
      
      event_time[treatment == 0] <- m$f0(runif(sum(treatment == 0)))
      event_time[treatment == 1] <- m$f1(runif(sum(treatment == 1)))
      
      data_sim <- data.frame(
        event_time = event_time,
        status = 1,
        treatment = treatment
      )
      
      fit <- cox_qr_estimate(data_sim, taus)
      beta_hat <- fit$coef_values$beta_hat
      
      boot <- bootstrap_cox_qr(data_sim, taus, B)
      
      CI_lower <- apply(boot, 2, quantile, probs = alpha/2)
      CI_upper <- apply(boot, 2, quantile, probs = 1 - alpha/2)
      
      covered <- (beta_true >= CI_lower) & (beta_true <= CI_upper)
      
      results[[idx]] <- data.frame(
        Model = m$name,
        sim = sim,
        tau = taus,
        beta_true = beta_true,
        beta_hat = beta_hat,
        covered = covered
      )
      
      idx <- idx + 1
    }
  }
  
  results <- do.call(rbind, results)
  
  coverage <- aggregate(covered ~ Model + tau, results, mean)
  
  write.csv(coverage, out_file, row.names = FALSE)
  
  invisible(coverage)
}


# -------------------------------------------------------------------------
# Coverage: Models VII-VIII
# -------------------------------------------------------------------------

run_coverage_VII_VIII <- function(n_sim = 2000,
                                  n = 300,
                                  p = 1/2,
                                  taus = taus,
                                  out_file = "simulations/coverage_models_VII_VIII.csv") {
  
  models <- get_models_VII_VIII(taus = taus, p = p)
  results <- vector("list", length(models) * n_sim)
  idx <- 1
  
  for (m in models) {
    
    aqe_true <- m$true_aqe
    
    for (sim in seq_len(n_sim)) {
      
      treatment <- rbinom(n, 1, p)
      X2 <- m$generate_X2(treatment)
      u <- runif(n)
      
      event_time <- m$generate_times(u, treatment, X2)
      
      data_sim <- data.frame(
        event_time = event_time,
        status = 1,
        treatment = treatment,
        X2 = X2
      )
      
      fit <- cox_qr_estimate_cov(data_sim, taus)
      aqe_hat <- fit$aqe_hat$AQE_hat
      
      boot <- bootstrap_cox_qr_cov(data_sim, taus, B)
      
      CI_lower <- apply(boot, 2, quantile, probs = alpha/2)
      CI_upper <- apply(boot, 2, quantile, probs = 1 - alpha/2)
      
      covered <- (aqe_true >= CI_lower) & (aqe_true <= CI_upper)
      
      results[[idx]] <- data.frame(
        Model = m$name,
        sim = sim,
        tau = taus,
        AQE_true = aqe_true,
        AQE_hat = aqe_hat,
        covered = covered
      )
      
      idx <- idx + 1
    }
  }
  
  results <- do.call(rbind, results)
  
  coverage <- aggregate(covered ~ Model + tau, results, mean)
  
  write.csv(coverage, out_file, row.names = FALSE)
  
  invisible(coverage)
}


# -------------------------------------------------------------------------
# Run
# -------------------------------------------------------------------------

run_coverage_I_VI(
  n_sim = 2000,
  n = 300,
  p = 1/3
)

run_coverage_VII_VIII(
  n_sim = 2000,
  n = 300,
  p = 1/2
)