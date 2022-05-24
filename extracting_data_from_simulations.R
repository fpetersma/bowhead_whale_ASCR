## ========================================================================== ##
## In this script, I extract all relevant info from the fits to the simulated ##
## data. These extracts are then stored in .csv files which are easier to     ##
## work with.                                                                 ##
## ========================================================================== ##

## =============================================================================
## 1. VARIABLE SOURCE LEVEL SIMULATIONS
## -----------------------------------------------------------------------------

## =============================================================================
## 1.1 FITS WITH VARIABLE SOURCE LEVEL AND MIXTURE BEARING MODEL
## -----------------------------------------------------------------------------

## Load data
load("Simulation output/fits_with_fixedSL=FALSE_USE_BEARINGS=2_1-100_on_simulations_with_fixedSL=FALSE.RData")
load("1000_simulations_fixedSL=FALSE.RData")

df_varSL <- data.frame(seed = 090596 + 1:100,
                       N_expected = NA,
                       N_realised = NA,
                       N_estimated = NA,
                       error_abs = NA,
                       error_rel = NA,
                       g0 = NA,
                       beta_r = NA,
                       sd_r = NA,
                       mu_s = NA,
                       sd_s = NA,
                       kappa_low = NA,
                       kappa_high = NA,
                       mix_bear = NA,
                       beta_intercept = NA,
                       beta_d2c = NA,
                       beta_d2c2 = NA)

## Loop through all fits and extract relevant info
for (i in 1:length(results)) {
  ## Extract the ith fit and simulation
  fit <- results[[i]]
  sim <- simulated_data[[i]]
  
  ## Extract the estimates
  
  df_varSL[i, "N_expected"] <- sim$N_expected
  df_varSL[i, "N_realised"] <- sim$N_realised
  df_varSL[i, "N_estimated"] <- fit$N_total
  df_varSL[i, "error_abs"] <- abs(df_varSL[i, "N_realised"] - 
                                     df_varSL[i, "N_estimated"])
  df_varSL[i, "error_rel"] <- df_varSL[i, "error_abs"] / df_varSL[i, "N_realised"]
  df_varSL[i, "g0"] <- fit$real_par["logit_g0"]
  df_varSL[i, "beta_r"] <- fit$real_par["log_beta_r"]
  df_varSL[i, "sd_r"] <- fit$real_par["log_sd_r"]
  df_varSL[i, "mu_s"] <- fit$real_par["log_mu_s"]
  df_varSL[i, "sd_s"] <- fit$real_par["log_sd_s"]
  df_varSL[i, "kappa_low"] <- fit$real_par["log_kappa_low"]
  df_varSL[i, "kappa_high"] <- fit$real_par["log_kappa_high"]
  df_varSL[i, "mix_bear"] <- fit$real_par["logit_mix_bear"]
  df_varSL[i, "beta_intercept"] <- fit$fit_res$par["(Intercept)"]
  df_varSL[i, "beta_d2c"] <- fit$fit_res$par["distance_to_coast"]
  df_varSL[i, "beta_d2c2"] <- fit$fit_res$par["distance_to_coast2"]
  # df_varSL[i, ""] <- fit$real_par[""]
}

## Create a summary data.frame
summary_varSL <- data.frame(quantity = colnames(df_varSL)[-1],
                            truth = c(NA, NA, NA, NA, NA,
                                      0.6,       # g0
                                      18,        # beta_r
                                      2.7,       # sd_r
                                      163,       # mu_s
                                      5,         # sd_s
                                      0.3,       # kappa_low
                                      37,        # kappa_high
                                      0.1,       # mix_bear
                                      -12,
                                      45,
                                      -53),
                            mean = apply(df_varSL[, -1], 2, mean),
                            sd = apply(df_varSL[, -1], 2, sd))
summary_varSL$cv <- summary_varSL$sd / abs(summary_varSL$mean) * 100
## -----------------------------------------------------------------------------

## =============================================================================
## 1.2 FITS WITH FIXED SOURCE LEVEL AND MIXTURE BEARING MODEL
## -----------------------------------------------------------------------------

## Load data
load("Simulation output/fits_with_fixedSL=TRUE_USE_BEARINGS=2_1-100_on_simulations_with_fixedSL=FALSE.RData")
load("1000_simulations_fixedSL=FALSE.RData")

df_varSL_fixedSL_fits <- data.frame(seed = 090596 + 1:100, 
                                    N_expected = NA,
                                    N_realised = NA,
                                    N_estimated = NA,
                                    error_abs = NA,
                                    error_rel = NA,
                                    g0 = NA,
                                    beta_r = NA,
                                    sd_r = NA,
                                    mu_s = NA,
                                    kappa_low = NA,
                                    kappa_high = NA,
                                    mix_bear = NA,
                                    beta_intercept = NA,
                                    beta_d2c = NA,
                                    beta_d2c2 = NA)

## Loop through all fits and extract relevant info
for (i in 1:length(results)) {
  ## Extract the ith fit and simulation
  fit <- results[[i]]
  sim <- simulated_data[[i]]
  
  ## Extract the estimates
  df_varSL_fixedSL_fits[i, "N_expected"] <- sim$N_expected
  df_varSL_fixedSL_fits[i, "N_realised"] <- sim$N_realised
  df_varSL_fixedSL_fits[i, "N_estimated"] <- fit$N_total
  df_varSL_fixedSL_fits[i, "error_abs"] <- abs(df_varSL_fixedSL_fits[i, "N_realised"] - 
                                                  df_varSL_fixedSL_fits[i, "N_estimated"])
  df_varSL_fixedSL_fits[i, "error_rel"] <- df_varSL_fixedSL_fits[i, "error_abs"] / df_varSL_fixedSL_fits[i, "N_realised"]
  df_varSL_fixedSL_fits[i, "g0"] <- fit$real_par["logit_g0"]
  df_varSL_fixedSL_fits[i, "beta_r"] <- fit$real_par["log_beta_r"]
  df_varSL_fixedSL_fits[i, "sd_r"] <- fit$real_par["log_sd_r"]
  df_varSL_fixedSL_fits[i, "mu_s"] <- fit$real_par["log_mu_s"]
  df_varSL_fixedSL_fits[i, "kappa_low"] <- fit$real_par["log_kappa_low"]
  df_varSL_fixedSL_fits[i, "kappa_high"] <- fit$real_par["log_kappa_high"]
  df_varSL_fixedSL_fits[i, "mix_bear"] <- fit$real_par["logit_mix_bear"]
  df_varSL_fixedSL_fits[i, "beta_intercept"] <- fit$fit_res$par["(Intercept)"]
  df_varSL_fixedSL_fits[i, "beta_d2c"] <- fit$fit_res$par["distance_to_coast"]
  df_varSL_fixedSL_fits[i, "beta_d2c2"] <- fit$fit_res$par["distance_to_coast2"]
  # df_varSL_fixedSL_fits[i, ""] <- fit$real_par[""]
}

## Create a summary data.frame
summary_varSL_fixedSL_fits <- data.frame(quantity = colnames(df_varSL_fixedSL_fits)[-1],
                                         truth = c(NA, NA, NA, NA, NA, 
                                                   0.6,       # g0
                                                   18,        # beta_r
                                                   2.7,       # sd_r
                                                   163,       # mu_s
                                                   0.3,       # kappa_low
                                                   37,        # kappa_high
                                                   0.1,       # mix_bear
                                                   -12,
                                                   45,
                                                   -53),
                                         mean = apply(df_varSL_fixedSL_fits[, -1], 2, mean),
                                         sd = apply(df_varSL_fixedSL_fits[, -1], 2, sd))
summary_varSL_fixedSL_fits$cv <- summary_varSL_fixedSL_fits$sd / 
  abs(summary_varSL_fixedSL_fits$mean) * 100
## -----------------------------------------------------------------------------

## =============================================================================
## 1.3 FITS WITH VARIABLE SOURCE LEVEL AND SIMPLE BEARING MODEL
## -----------------------------------------------------------------------------

## Load data
load("Simulation output/fits_with_fixedSL=FALSE_USE_BEARINGS=1_1-100_on_simulations_with_fixedSL=FALSE.RData")
load("1000_simulations_fixedSL=FALSE.RData")

df_varSL_simple_bearing <- data.frame(seed = 090596 + 1:100,
                                      N_expected = NA,
                                      N_realised = NA,
                                      N_estimated = NA,
                                      error_abs = NA,
                                      error_rel = NA,
                                      g0 = NA,
                                      beta_r = NA,
                                      sd_r = NA,
                                      mu_s = NA,
                                      sd_s = NA,
                                      kappa = NA,
                                      beta_intercept = NA,
                                      beta_d2c = NA,
                                      beta_d2c2 = NA)

## Loop through all fits and extract relevant info
for (i in 1:length(results)) {
  ## Extract the ith fit and simulation
  fit <- results[[i]]
  sim <- simulated_data[[i]]
  
  ## Extract the estimates
  df_varSL_simple_bearing[i, "N_expected"] <- sim$N_expected
  df_varSL_simple_bearing[i, "N_realised"] <- sim$N_realised
  df_varSL_simple_bearing[i, "N_estimated"] <- fit$N_total
  df_varSL_simple_bearing[i, "error_abs"] <- abs(df_varSL_simple_bearing[i, "N_realised"] - 
                                                    df_varSL_simple_bearing[i, "N_estimated"])
  df_varSL_simple_bearing[i, "error_rel"] <- df_varSL_simple_bearing[i, "error_abs"] / 
    df_varSL_simple_bearing[i, "N_realised"]
  df_varSL_simple_bearing[i, "g0"] <- fit$real_par["logit_g0"]
  df_varSL_simple_bearing[i, "beta_r"] <- fit$real_par["log_beta_r"]
  df_varSL_simple_bearing[i, "sd_r"] <- fit$real_par["log_sd_r"]
  df_varSL_simple_bearing[i, "mu_s"] <- fit$real_par["log_mu_s"]
  df_varSL_simple_bearing[i, "sd_s"] <- fit$real_par["log_sd_s"]
  df_varSL_simple_bearing[i, "kappa"] <- fit$real_par["log_kappa"]
  df_varSL_simple_bearing[i, "beta_intercept"] <- fit$fit_res$par["(Intercept)"]
  df_varSL_simple_bearing[i, "beta_d2c"] <- fit$fit_res$par["distance_to_coast"]
  df_varSL_simple_bearing[i, "beta_d2c2"] <- fit$fit_res$par["distance_to_coast2"]
}

## Create a summary data.frame
summary_varSL_simple_bearing <- data.frame(quantity = colnames(df_varSL_simple_bearing)[-1],
                                           truth = c(NA, NA, NA, NA, NA, 
                                                     0.6,       # g0
                                                     18,        # beta_r
                                                     2.7,       # sd_r
                                                     163,       # mu_s
                                                     5,         # sd_s
                                                     NA,       # kappa
                                                     -12,
                                                     45,
                                                     -53),
                                           mean = apply(df_varSL_simple_bearing[, -1], 2, mean),
                                           sd = apply(df_varSL_simple_bearing[, -1], 2, sd))
summary_varSL_simple_bearing$cv <- summary_varSL_simple_bearing$sd / 
  abs(summary_varSL_simple_bearing$mean) * 100

## =============================================================================
## 1.4 FITS WITH VARIABLE SOURCE LEVEL AND WITHOUT BEARINGS
## -----------------------------------------------------------------------------

## Load data
load("Simulation output/fits_with_fixedSL=FALSE_USE_BEARINGS=0_1-100_on_simulations_with_fixedSL=FALSE.RData")
load("1000_simulations_fixedSL=FALSE.RData")

df_varSL_no_bearings <- data.frame(seed = 090596 + 1:100,
                                   N_expected = NA,
                                   N_realised = NA,
                                   N_estimated = NA,
                                   error_abs = NA,
                                   error_rel = NA,
                                   g0 = NA,
                                   beta_r = NA,
                                   sd_r = NA,
                                   mu_s = NA,
                                   sd_s = NA,
                                   beta_intercept = NA,
                                   beta_d2c = NA,
                                   beta_d2c2 = NA)

## Loop through all fits and extract relevant info
for (i in 1:length(results)) {
  ## Extract the ith fit and simulation
  fit <- results[[i]]
  sim <- simulated_data[[i]]
  
  ## Extract the estimates
  df_varSL_no_bearings[i, "N_expected"] <- sim$N_expected
  df_varSL_no_bearings[i, "N_realised"] <- sim$N_realised
  df_varSL_no_bearings[i, "N_estimated"] <- fit$N_total
  df_varSL_no_bearings[i, "error_abs"] <- abs(df_varSL_no_bearings[i, "N_realised"] - 
                                                 df_varSL_no_bearings[i, "N_estimated"])
  df_varSL_no_bearings[i, "error_rel"] <- df_varSL_no_bearings[i, "error_abs"] / 
    df_varSL_no_bearings[i, "N_realised"]
  df_varSL_no_bearings[i, "g0"] <- fit$real_par["logit_g0"]
  df_varSL_no_bearings[i, "beta_r"] <- fit$real_par["log_beta_r"]
  df_varSL_no_bearings[i, "sd_r"] <- fit$real_par["log_sd_r"]
  df_varSL_no_bearings[i, "mu_s"] <- fit$real_par["log_mu_s"]
  df_varSL_no_bearings[i, "sd_s"] <- fit$real_par["log_sd_s"]
  df_varSL_no_bearings[i, "beta_intercept"] <- fit$fit_res$par["(Intercept)"]
  df_varSL_no_bearings[i, "beta_d2c"] <- fit$fit_res$par["distance_to_coast"]
  df_varSL_no_bearings[i, "beta_d2c2"] <- fit$fit_res$par["distance_to_coast2"]
}

## Create a summary data.frame
summary_varSL_no_bearings <- data.frame(quantity = colnames(df_varSL_no_bearings)[-1],
                                        truth = c(NA, NA, NA, NA, NA, 
                                                  0.6,       # g0
                                                  18,        # beta_r
                                                  2.7,       # sd_r
                                                  163,       # mu_s
                                                  5,         # sd_s
                                                  -12,
                                                  45,
                                                  -53),
                                        mean = apply(df_varSL_no_bearings[, -1], 2, mean),
                                        sd = apply(df_varSL_no_bearings[, -1], 2, sd))
summary_varSL_no_bearings$cv <- summary_varSL_no_bearings$sd / 
  abs(summary_varSL_no_bearings$mean) * 100

## =============================================================================
## 1.5 FITS WITH VARIABLE SOURCE LEVEL, MIXTURE BEARING MODEL AND CONSTANT D
## -----------------------------------------------------------------------------

## Load data
load("Simulation output/fits_with_fixedSL=FALSE_USE_BEARINGS=2_constant_density_1-100_on_simulations_with_fixedSL=FALSE.RData")
load("1000_simulations_fixedSL=FALSE.RData")

df_varSL_constant_D <- data.frame(seed = 090596 + 1:100,
                                  N_expected = NA,
                                  N_realised = NA,
                                  N_estimated = NA,
                                  error_abs = NA,
                                  error_rel = NA,
                                  g0 = NA,
                                  beta_r = NA,
                                  sd_r = NA,
                                  mu_s = NA,
                                  sd_s = NA,
                                  kappa_low = NA,
                                  kappa_high = NA,
                                  mix_bear = NA,
                                  beta_intercept = NA)

## Loop through all fits and extract relevant info
for (i in 1:length(results)) {
  ## Extract the ith fit and simulation
  fit <- results[[i]]
  sim <- simulated_data[[i]]
  
  ## Extract the estimates
  df_varSL_constant_D[i, "N_expected"] <- sim$N_expected
  df_varSL_constant_D[i, "N_realised"] <- sim$N_realised
  df_varSL_constant_D[i, "N_estimated"] <- fit$N_total
  df_varSL_constant_D[i, "error_abs"] <- abs(df_varSL_constant_D[i, "N_realised"] - 
                                                df_varSL_constant_D[i, "N_estimated"])
  df_varSL_constant_D[i, "error_rel"] <- df_varSL_constant_D[i, "error_abs"] / 
    df_varSL_constant_D[i, "N_realised"]
  df_varSL_constant_D[i, "g0"] <- fit$real_par["logit_g0"]
  df_varSL_constant_D[i, "beta_r"] <- fit$real_par["log_beta_r"]
  df_varSL_constant_D[i, "sd_r"] <- fit$real_par["log_sd_r"]
  df_varSL_constant_D[i, "mu_s"] <- fit$real_par["log_mu_s"]
  df_varSL_constant_D[i, "sd_s"] <- fit$real_par["log_sd_s"]
  df_varSL_constant_D[i, "kappa_low"] <- fit$real_par["log_kappa_low"]
  df_varSL_constant_D[i, "kappa_high"] <- fit$real_par["log_kappa_high"]
  df_varSL_constant_D[i, "mix_bear"] <- fit$real_par["logit_mix_bear"]
  df_varSL_constant_D[i, "beta_intercept"] <- fit$fit_res$par["(Intercept)"]
}

## Create a summary data.frame
summary_varSL_constant_D <- data.frame(quantity = colnames(df_varSL_constant_D)[-1],
                                       truth = c(NA, NA, NA, NA, NA, 
                                                 0.6,       # g0
                                                 18,        # beta_r
                                                 2.7,       # sd_r
                                                 163,       # mu_s
                                                 5,         # sd_s
                                                 0.3,       # kappa_low
                                                 37,        # kappa_high
                                                 0.1,       # mix_bear
                                                 -NA),      # beta_intercept
                                       mean = apply(df_varSL_constant_D[, -1], 2, mean),
                                       sd = apply(df_varSL_constant_D[, -1], 2, sd))
summary_varSL_constant_D$cv <- summary_varSL_constant_D$sd / 
  abs(summary_varSL_constant_D$mean)  * 100

## =============================================================================
## 2 FIXED SOURCE LEVEL SIMULATIONS
## -----------------------------------------------------------------------------

## =============================================================================
## 2.1 FITS WITH FIXED SOURCE LEVEL AND MIXTURE BEARING MODEL
## -----------------------------------------------------------------------------

## Load data
load("Simulation output/fits_with_fixedSL=TRUE_USE_BEARINGS=2_1-100_on_simulations_with_fixedSL=TRUE.RData")
load("1000_simulations_fixedSL=TRUE.RData")

df_fixedSL <- data.frame(seed = 090596 + 1:100, 
                         N_expected = NA,
                         N_realised = NA,
                         N_estimated = NA,
                         error_abs = NA,
                         error_rel = NA,
                         g0 = NA,
                         beta_r = NA,
                         sd_r = NA,
                         mu_s = NA,
                         kappa_low = NA,
                         kappa_high = NA,
                         mix_bear = NA,
                         beta_intercept = NA,
                         beta_d2c = NA,
                         beta_d2c2 = NA)

## Loop through all fits and extract relevant info
for (i in 1:length(results)) {
  ## Extract the ith fit and simulation
  fit <- results[[i]]
  sim <- simulated_data[[i]]
  
  ## Extract the estimates
  df_fixedSL[i, "N_expected"] <- sim$N_expected
  df_fixedSL[i, "N_realised"] <- sim$N_realised
  df_fixedSL[i, "N_estimated"] <- fit$N_total
  df_fixedSL[i, "error_abs"] <- abs(df_fixedSL[i, "N_realised"] - 
                                       df_fixedSL[i, "N_estimated"])
  df_fixedSL[i, "error_rel"] <- df_fixedSL[i, "error_abs"] / df_fixedSL[i, "N_realised"]
  df_fixedSL[i, "g0"] <- fit$real_par["logit_g0"]
  df_fixedSL[i, "beta_r"] <- fit$real_par["log_beta_r"]
  df_fixedSL[i, "sd_r"] <- fit$real_par["log_sd_r"]
  df_fixedSL[i, "mu_s"] <- fit$real_par["log_mu_s"]
  df_fixedSL[i, "kappa_low"] <- fit$real_par["log_kappa_low"]
  df_fixedSL[i, "kappa_high"] <- fit$real_par["log_kappa_high"]
  df_fixedSL[i, "mix_bear"] <- fit$real_par["logit_mix_bear"]
  df_fixedSL[i, "beta_intercept"] <- fit$fit_res$par["(Intercept)"]
  df_fixedSL[i, "beta_d2c"] <- fit$fit_res$par["distance_to_coast"]
  df_fixedSL[i, "beta_d2c2"] <- fit$fit_res$par["distance_to_coast2"]
  # df_fixedSL[i, ""] <- fit$real_par[""]
}

## Create a summary data.frame
summary_fixedSL <- data.frame(quantity = colnames(df_fixedSL)[-1],
                              truth = c(NA, NA, NA, NA, NA, 
                                        0.6,       # g0
                                        14.5,      # beta_r
                                        4.5,       # sd_r
                                        155,       # mu_s
                                        0.3,       # kappa_low
                                        37,        # kappa_high
                                        0.1,       # mix_bear
                                        -16,
                                        57,
                                        -68.5),
                              mean = apply(df_fixedSL[, -1], 2, mean),
                              sd = apply(df_fixedSL[, -1], 2, sd))
summary_fixedSL$cv <- summary_fixedSL$sd / abs(summary_fixedSL$mean)  * 100
## -----------------------------------------------------------------------------

## =============================================================================
## 2.2 FITS WITH VARIABLE SOURCE LEVEL AND MIXTURE BEARING MODEL
## -----------------------------------------------------------------------------

## Load data
load("Simulation output/fits_with_fixedSL=FALSE_USE_BEARINGS=2_1-100_on_simulations_with_fixedSL=TRUE_density_pars=c(-16,57,-68.5).RData")
load("1000_simulations_fixedSL=TRUE.RData")

df_fixedSL_varSL_fits <- data.frame(seed = 090596 + 1:100, 
                                    N_expected = NA,
                                    N_realised = NA,
                                    N_estimated = NA,
                                    error_abs = NA,
                                    error_rel = NA,
                                    g0 = NA,
                                    beta_r = NA,
                                    sd_r = NA,
                                    mu_s = NA,
                                    sd_s = NA,
                                    kappa_low = NA,
                                    kappa_high = NA,
                                    mix_bear = NA,
                                    beta_intercept = NA,
                                    beta_d2c = NA,
                                    beta_d2c2 = NA)

## Loop through all fits and extract relevant info
for (i in 1:length(results)) {
  ## Extract the ith fit and simulation
  fit <- results[[i]]
  sim <- simulated_data[[i]]
  
  ## Extract the estimates
  df_fixedSL_varSL_fits[i, "N_expected"] <- sim$N_expected
  df_fixedSL_varSL_fits[i, "N_realised"] <- sim$N_realised
  df_fixedSL_varSL_fits[i, "N_estimated"] <- fit$N_total
  df_fixedSL_varSL_fits[i, "error_abs"] <- abs(df_fixedSL_varSL_fits[i, "N_realised"] - 
                                                  df_fixedSL_varSL_fits[i, "N_estimated"])
  df_fixedSL_varSL_fits[i, "error_rel"] <- df_fixedSL_varSL_fits[i, "error_abs"] / 
    df_fixedSL_varSL_fits[i, "N_realised"]
  df_fixedSL_varSL_fits[i, "g0"] <- fit$real_par["logit_g0"]
  df_fixedSL_varSL_fits[i, "beta_r"] <- fit$real_par["log_beta_r"]
  df_fixedSL_varSL_fits[i, "sd_r"] <- fit$real_par["log_sd_r"]
  df_fixedSL_varSL_fits[i, "mu_s"] <- fit$real_par["log_mu_s"]
  df_fixedSL_varSL_fits[i, "sd_s"] <- fit$real_par["log_sd_s"]
  df_fixedSL_varSL_fits[i, "kappa_low"] <- fit$real_par["log_kappa_low"]
  df_fixedSL_varSL_fits[i, "kappa_high"] <- fit$real_par["log_kappa_high"]
  df_fixedSL_varSL_fits[i, "mix_bear"] <- fit$real_par["logit_mix_bear"]
  df_fixedSL_varSL_fits[i, "beta_intercept"] <- fit$fit_res$par["(Intercept)"]
  df_fixedSL_varSL_fits[i, "beta_d2c"] <- fit$fit_res$par["distance_to_coast"]
  df_fixedSL_varSL_fits[i, "beta_d2c2"] <- fit$fit_res$par["distance_to_coast2"]
  # df_fixedSL[i, ""] <- fit$real_par[""]
}

## Create a summary data.frame
summary_fixedSL_varSL_fits <- data.frame(quantity = colnames(df_fixedSL_varSL_fits)[-1],
                                         truth = c(NA, NA, NA, NA, NA, 
                                                   0.6,       # g0
                                                   14.5,      # beta_r
                                                   4.5,       # sd_r
                                                   155,       # mu_s
                                                   NA,        # sd_s
                                                   0.3,       # kappa_low
                                                   37,        # kappa_high
                                                   0.1,       # mix_bear
                                                   -16,
                                                   57,
                                                   -68.5),
                                         mean = apply(df_fixedSL_varSL_fits[, -1], 2, mean),
                                         sd = apply(df_fixedSL_varSL_fits[, -1], 2, sd))
summary_fixedSL_varSL_fits$cv <- summary_fixedSL_varSL_fits$sd / 
  abs(summary_fixedSL_varSL_fits$mean) * 100
## -----------------------------------------------------------------------------

## =============================================================================
## 2.3 FITS WITH FIXED SOURCE LEVEL AND SIMPLE BEARING MODEL
## -----------------------------------------------------------------------------

## Load data
load("Simulation output/fits_with_fixedSL=TRUE_USE_BEARINGS=1_1-100_on_simulations_with_fixedSL=TRUE.RData")
load("1000_simulations_fixedSL=TRUE.RData")

df_fixedSL_simple_bearing <- data.frame(seed = 090596 + 1:100, 
                                        N_expected = NA,
                                        N_realised = NA,
                                        N_estimated = NA,
                                        error_abs = NA,
                                        error_rel = NA,
                                        g0 = NA,
                                        beta_r = NA,
                                        sd_r = NA,
                                        mu_s = NA,
                                        kappa = NA,
                                        beta_intercept = NA,
                                        beta_d2c = NA,
                                        beta_d2c2 = NA)

## Loop through all fits and extract relevant info
for (i in 1:length(results)) {
  ## Extract the ith fit and simulation
  fit <- results[[i]]
  sim <- simulated_data[[i]]
  
  ## Extract the estimates
  df_fixedSL_simple_bearing[i, "N_expected"] <- sim$N_expected
  df_fixedSL_simple_bearing[i, "N_realised"] <- sim$N_realised
  df_fixedSL_simple_bearing[i, "N_estimated"] <- fit$N_total
  df_fixedSL_simple_bearing[i, "error_abs"] <- abs(df_fixedSL_simple_bearing[i, "N_realised"] - 
                                                      df_fixedSL_simple_bearing[i, "N_estimated"])
  df_fixedSL_simple_bearing[i, "error_rel"] <- df_fixedSL_simple_bearing[i, "error_abs"] / 
    df_fixedSL_simple_bearing[i, "N_realised"]
  df_fixedSL_simple_bearing[i, "g0"] <- fit$real_par["logit_g0"]
  df_fixedSL_simple_bearing[i, "beta_r"] <- fit$real_par["log_beta_r"]
  df_fixedSL_simple_bearing[i, "sd_r"] <- fit$real_par["log_sd_r"]
  df_fixedSL_simple_bearing[i, "mu_s"] <- fit$real_par["log_mu_s"]
  df_fixedSL_simple_bearing[i, "kappa"] <- fit$real_par["log_kappa"]
  df_fixedSL_simple_bearing[i, "beta_intercept"] <- fit$fit_res$par["(Intercept)"]
  df_fixedSL_simple_bearing[i, "beta_d2c"] <- fit$fit_res$par["distance_to_coast"]
  df_fixedSL_simple_bearing[i, "beta_d2c2"] <- fit$fit_res$par["distance_to_coast2"]
}

## Create a summary data.frame
summary_fixedSL_simple_bearing <- data.frame(quantity = colnames(df_fixedSL_simple_bearing)[-1],
                                             truth = c(NA, NA, NA, NA, NA, 
                                                       0.6,       # g0
                                                       14.5,      # beta_r
                                                       4.5,       # sd_r
                                                       155,       # mu_s
                                                       NA,        # kappa
                                                       -16,
                                                       57,
                                                       -68.5),
                                             mean = apply(df_fixedSL_simple_bearing[, -1], 2, mean),
                                             sd = apply(df_fixedSL_simple_bearing[, -1], 2, sd))
summary_fixedSL_simple_bearing$cv <- summary_fixedSL_simple_bearing$sd / 
  abs(summary_fixedSL_simple_bearing$mean) * 100
## -----------------------------------------------------------------------------

## =============================================================================
## 2.4 FITS WITH FIXED SOURCE LEVEL AND WITHOUT BEARINGS
## -----------------------------------------------------------------------------

## Load data
load("Simulation output/fits_with_fixedSL=TRUE_USE_BEARINGS=0_1-100_on_simulations_with_fixedSL=TRUE.RData")
load("1000_simulations_fixedSL=TRUE.RData")

df_fixedSL_no_bearings <- data.frame(seed = 090596 + 1:100, 
                                     N_expected = NA,
                                     N_realised = NA,
                                     N_estimated = NA,
                                     error_abs = NA,
                                     error_rel = NA,
                                     g0 = NA,
                                     beta_r = NA,
                                     sd_r = NA,
                                     mu_s = NA,
                                     beta_intercept = NA,
                                     beta_d2c = NA,
                                     beta_d2c2 = NA)

## Loop through all fits and extract relevant info
for (i in 1:length(results)) {
  ## Extract the ith fit and simulation
  fit <- results[[i]]
  sim <- simulated_data[[i]]
  
  ## Extract the estimates
  df_fixedSL_no_bearings[i, "N_expected"] <- sim$N_expected
  df_fixedSL_no_bearings[i, "N_realised"] <- sim$N_realised
  df_fixedSL_no_bearings[i, "N_estimated"] <- fit$N_total
  df_fixedSL_no_bearings[i, "error_abs"] <- abs(df_fixedSL_no_bearings[i, "N_realised"] - 
                                                   df_fixedSL_no_bearings[i, "N_estimated"])
  df_fixedSL_no_bearings[i, "error_rel"] <- df_fixedSL_no_bearings[i, "error_abs"] / 
    df_fixedSL_no_bearings[i, "N_realised"]
  df_fixedSL_no_bearings[i, "g0"] <- fit$real_par["logit_g0"]
  df_fixedSL_no_bearings[i, "beta_r"] <- fit$real_par["log_beta_r"]
  df_fixedSL_no_bearings[i, "sd_r"] <- fit$real_par["log_sd_r"]
  df_fixedSL_no_bearings[i, "mu_s"] <- fit$real_par["log_mu_s"]
  df_fixedSL_no_bearings[i, "beta_intercept"] <- fit$fit_res$par["(Intercept)"]
  df_fixedSL_no_bearings[i, "beta_d2c"] <- fit$fit_res$par["distance_to_coast"]
  df_fixedSL_no_bearings[i, "beta_d2c2"] <- fit$fit_res$par["distance_to_coast2"]
}

## Create a summary data.frame
summary_fixedSL_no_bearings <- data.frame(quantity = colnames(df_fixedSL_no_bearings)[-1],
                                          truth = c(NA, NA, NA, NA, NA, 
                                                    0.6,       # g0
                                                    14.5,      # beta_r
                                                    4.5,       # sd_r
                                                    155,       # mu_s
                                                    -16,
                                                    57,
                                                    -68.5),
                                          mean = apply(df_fixedSL_no_bearings[, -1], 2, mean),
                                          sd = apply(df_fixedSL_no_bearings[, -1], 2, sd))
summary_fixedSL_no_bearings$cv <- summary_fixedSL_no_bearings$sd / 
  abs(summary_fixedSL_no_bearings$mean) * 100
## -----------------------------------------------------------------------------

## =============================================================================
## 2.5 FITS WITH FIXED SOURCE LEVEL, MIXTURE BEARING MODEL AND CONSTANT D
## -----------------------------------------------------------------------------

## Load data
load("Simulation output/fits_with_fixedSL=TRUE_USE_BEARINGS=2_constant_density_1-100_on_simulations_with_fixedSL=TRUE.RData")
load("1000_simulations_fixedSL=TRUE.RData")

df_fixedSL_constant_D <- data.frame(seed = 090596 + 1:100, 
                                    N_expected = NA,
                                    N_realised = NA,
                                    N_estimated = NA,
                                    error_abs = NA,
                                    error_rel = NA,
                                    g0 = NA,
                                    beta_r = NA,
                                    sd_r = NA,
                                    mu_s = NA,
                                    kappa_low = NA,
                                    kappa_high = NA,
                                    mix_bear = NA,
                                    beta_intercept = NA)

## Loop through all fits and extract relevant info
for (i in 1:length(results)) {
  ## Extract the ith fit and simulation
  fit <- results[[i]]
  sim <- simulated_data[[i]]
  
  ## Extract the estimates
  df_fixedSL_constant_D[i, "N_expected"] <- sim$N_expected
  df_fixedSL_constant_D[i, "N_realised"] <- sim$N_realised
  df_fixedSL_constant_D[i, "N_estimated"] <- fit$N_total
  df_fixedSL_constant_D[i, "error_abs"] <- abs(df_fixedSL_constant_D[i, "N_realised"] - 
                                                  df_fixedSL_constant_D[i, "N_estimated"])
  df_fixedSL_constant_D[i, "error_rel"] <- df_fixedSL_constant_D[i, "error_abs"] / 
    df_fixedSL_constant_D[i, "N_realised"]
  df_fixedSL_constant_D[i, "g0"] <- fit$real_par["logit_g0"]
  df_fixedSL_constant_D[i, "beta_r"] <- fit$real_par["log_beta_r"]
  df_fixedSL_constant_D[i, "sd_r"] <- fit$real_par["log_sd_r"]
  df_fixedSL_constant_D[i, "mu_s"] <- fit$real_par["log_mu_s"]
  df_fixedSL_constant_D[i, "kappa_low"] <- fit$real_par["log_kappa_low"]
  df_fixedSL_constant_D[i, "kappa_high"] <- fit$real_par["log_kappa_high"]
  df_fixedSL_constant_D[i, "mix_bear"] <- fit$real_par["logit_mix_bear"]
  df_fixedSL_constant_D[i, "beta_intercept"] <- fit$fit_res$par["(Intercept)"]
  # df_fixedSL_constant_D[i, ""] <- fit$real_par[""]
}

## Create a summary data.frame
summary_fixedSL_constant_D <- data.frame(quantity = colnames(df_fixedSL_constant_D)[-1],
                                         truth = c(NA, NA, NA, NA, NA, 
                                                   0.6,       # g0
                                                   14.5,      # beta_r
                                                   4.5,       # sd_r
                                                   155,       # mu_s
                                                   0.3,       # kappa_low
                                                   37,        # kappa_high
                                                   0.1,       # mix_bear
                                                   NA),       # beta_intercept
                                         mean = apply(df_fixedSL_constant_D[, -1], 2, mean),
                                         sd = apply(df_fixedSL_constant_D[, -1], 2, sd))
summary_fixedSL_constant_D$cv <- summary_fixedSL_constant_D$sd / 
  abs(summary_fixedSL_constant_D$mean) * 100
## -----------------------------------------------------------------------------

## Save df's as an RData file for use in other documents
save(list = c("df_fixedSL",
              "df_fixedSL_constant_D", 
              "df_fixedSL_no_bearings", 
              "df_fixedSL_varSL_fits", 
              "df_fixedSL_simple_bearing", 
              "df_varSL", 
              "df_varSL_constant_D", 
              "df_varSL_no_bearings", 
              "df_varSL_fixedSL_fits", 
              "df_varSL_simple_bearing"),
     file = "Simulation output/simulation_results_full.RData")

## Save summaries as an RData file for use in other documents
save(list = c("summary_fixedSL",
              "summary_fixedSL_constant_D", 
              "summary_fixedSL_no_bearings", 
              "summary_fixedSL_varSL_fits", 
              "summary_fixedSL_simple_bearing", 
              "summary_varSL", 
              "summary_varSL_constant_D", 
              "summary_varSL_no_bearings", 
              "summary_varSL_fixedSL_fits", 
              "summary_varSL_simple_bearing"),
     file = "Simulation output/simulation_results_summaries.RData")

## =============================================================================
## Get data ready for box plots and make plots
## -----------------------------------------------------------------------------

N_full <- data.frame(N_true_con = df_full$true_N,
                     N_variable = df_full$estimated_N,
                     N_constant = df_varSL_fixedSL_fits$estimated_N)
N_simple <- data.frame(N_true_var = df_fixedSL$true_N,
                       N_est_cons = df_fixedSL$estimated_N,
                       N_est_kappa = df_single_kappa$estimated_N)
N_full <- reshape2::melt(N_full)
N_simple <- reshape2::melt(N_simple)
library(ggplot2)
g_full <- ggplot(N_full) +
  geom_boxplot(aes(x = variable, y = value), fill = "lightgrey") +
  theme_bw() + ylab(element_blank()) + xlab(element_blank()) +
  scale_x_discrete(labels= c("Variable D\nsimulations", 
                             "Variable D\nfits",
                             "Constant D\nfits"))

g_simple <- ggplot(N_simple) +
  geom_boxplot(aes(x = variable, y = value), fill = "lightgrey") +
  theme_bw() + ylab(element_blank()) + xlab(element_blank()) +
  scale_x_discrete(labels= c("Constant D\nsimulations", 
                             "Constant D\nfits",
                             "Single bearing\nerror fits")) 



library(gridExtra)
grid.arrange(g_full, g_simple, nrow = 1, ncol = 2, left = "Estimated N", 
             bottom = "Model")


