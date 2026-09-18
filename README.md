# Cox-based quantile regression estimator

R code accompanying the manuscript:

Fernández-Iglesias et al.  
"A semiparametric estimator for quantile regression based on Cox proportional hazards models"

The repository contains the implementation of the proposed Cox-based quantile regression estimator and the code used for the simulation studies presented in the manuscript and Supplementary Material.

## Files

**cox_qr_estimator.R**  
Functions implementing the proposed Cox-based quantile regression estimator. The estimator allows for right-censored observations through the event indicator.

**simulation_models.R**  
Data-generating mechanisms used in the simulation study.

**simulation_errors.R**  
Code for the Monte Carlo simulation study evaluating estimation error under the different simulation models and sample-size configurations considered in the manuscript.

**simulation_coverage.R**  
Code for the simulation study evaluating the coverage of bootstrap confidence intervals.

**simulation_censoring.R**  
Code for the additional simulation study with 20% and 40% right censoring. Cox-based quantile regression is compared with censored quantile regression using the Peng-Huang and Portnoy methods.

## Requirements

The analyses were performed in R. The following packages are required:

- `survival`
- `quantreg`

Additional packages may be required for producing the figures reported in the manuscript.

## Simulation settings

The main simulation study considers the quantile levels 0.25, 0.50, and 0.75 and the following sample-size configurations:

- n = 100, p = 0.5
- n = 300, p = 0.5
- n = 300, p = 1/3

The results reported in the manuscript are based on 2,000 Monte Carlo replications. Bootstrap confidence intervals are based on 2,000 bootstrap samples.

The additional censoring simulations consider censoring proportions of 20% and 40%, in addition to the corresponding uncensored scenarios.

## Running the code

The scripts assume that all R files are located in the same directory. The simulation scripts source `cox_qr_estimator.R` and `simulation_models.R` as required.

The full simulation study is computationally intensive. The number of Monte Carlo or bootstrap replications can be reduced in the corresponding scripts for testing purposes.
