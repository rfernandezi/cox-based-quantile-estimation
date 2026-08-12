# -------------------------------------------------------------------------
# Monte Carlo simulations for Cox-based semiparametric quantile regression
# -------------------------------------------------------------------------
# Produces CSV files with estimation errors:
#   - Models I-VI: beta(tau) = Q_tau(Y | treatment = 1) - Q_tau(Y | treatment = 0)
#   - Models VII-VIII: AQE(tau) = E[beta(tau, X2)]
# -------------------------------------------------------------------------

library(survival)

source("R/cox_qr_estimator.R")
source("R/simulation_models.R")

set.seed(123)

taus <- c(0.25, 0.50, 0.75)

save_csv <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write.csv(x, path, row.names = FALSE)
}


# -------------------------------------------------------------------------
# Models I-VI: no covariates
# -------------------------------------------------------------------------

run_models_I_VI <- function(n_sim = 2000,
                            n = 300,
                            p = 1 / 3,
                            taus = c(0.25, 0.50, 0.75),
                            out_file = "simulations/errors_models_I_VI.csv") {
  
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
      
      beta_hat <- cox_qr_estimate(data_sim, taus)$coef_values$beta_hat
      
      results[[idx]] <- data.frame(
        Model = m$name,
        sim = sim,
        tau = taus,
        beta_true = beta_true,
        beta_hat = beta_hat,
        error = beta_hat - beta_true
      )
      
      idx <- idx + 1
    }
  }
  
  results <- do.call(rbind, results)
  save_csv(results, out_file)
  invisible(results)
}


# -------------------------------------------------------------------------
# Models VII-VIII: one covariate, AQE target
# -------------------------------------------------------------------------

run_models_VII_VIII <- function(n_sim = 2000,
                                n = 300,
                                p = 1 / 2,
                                taus = c(0.25, 0.50, 0.75),
                                out_file = "simulations/errors_models_VII_VIII.csv") {
  
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
      
      results[[idx]] <- data.frame(
        Model = m$name,
        sim = sim,
        tau = taus,
        AQE_true = aqe_true,
        AQE_hat = aqe_hat,
        error = aqe_hat - aqe_true
      )
      
      idx <- idx + 1
    }
  }
  
  results <- do.call(rbind, results)
  save_csv(results, out_file)
  invisible(results)
}


# -------------------------------------------------------------------------
# Run simulations
# -------------------------------------------------------------------------

FAST <- TRUE

if (FAST) {
  n_sim_I_VI <- 200
  n_sim_VII_VIII <- 200
} else {
  n_sim_I_VI <- 2000
  n_sim_VII_VIII <- 2000
}

errors_I_VI <- run_models_I_VI(
  n_sim = n_sim_I_VI,
  n = 300,
  p = 1 / 3,
  taus = taus,
  out_file = "simulations/errors_models_I_VI.csv"
)

errors_VII_VIII <- run_models_VII_VIII(
  n_sim = n_sim_VII_VIII,
  n = 300,
  p = 1 / 2,
  taus = taus,
  out_file = "simulations/errors_models_VII_VIII.csv"
)