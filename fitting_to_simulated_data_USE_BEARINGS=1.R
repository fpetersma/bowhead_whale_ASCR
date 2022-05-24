## ========================================================================== ##
## This script is run to fit data to the simulated                            ##
## ========================================================================== ##

## =============================================================================
## 1. LOAD THE SIMULATED DATA
## -----------------------------------------------------------------------------
load("1000_simulations_fixedSL=TRUE.RData")
## -----------------------------------------------------------------------------

## =============================================================================
## 2. LOOP THROUGH THE DATA, PREPARE AND FIT
## -----------------------------------------------------------------------------

## Number of data sets to include
n <- 100

## Use apply for more efficient looping
results <- lapply(simulated_data[1:2], function(dat) {
  ## Below is copied from mainRcpp.R
  
  ## Libraries
  library(Rcpp)       # for combining R and cpp
  library(parallel)   # for parallel processing
  library(mgcv)
  
  ## Constants
  A_s <- 3
  FIXED_SL <- FALSE
  USE_BEARINGS <- 1
  
  ## ===========================================================================
  ## 2.1 SET UP DENSITY MODEL
  ## ---------------------------------------------------------------------------
  ## Fit using the correct density model
  f_density <- D ~ distance_to_coast + distance_to_coast2
  
  ## Create GAM object with scaled covariates
  gam_fit <- gam(formula = f_density,
                 data = data.frame(D = 0, dat$covariates))
  
  ## Create data object
  data_obj <- list(Y_rec = dat$bearings_rad,      # in radians
                   Y_grid = dat$grid_bearings,    # in radians
                   X = dat$distances,
                   W = dat$det_hist,
                   R = dat$received_levels,
                   A_x = dat$A_x,
                   f_density = f_density,
                   trunc_level = dat$trunc_level,
                   gam_fit = gam_fit,
                   n_det = ncol(dat$distances),
                   n_grid = nrow(dat$distances),
                   n_call = nrow(dat$det_hist),
                   FIXED_SL = FIXED_SL,
                   USE_BEARINGS = USE_BEARINGS) # create data
  
  if (!FIXED_SL) {
    data_obj$A_s <- A_s
    data_obj$S <- seq(from = 100 + A_s / 2, to = 220 - A_s / 2, by = A_s) # for midpoint rule
    data_obj$n_sl <- length(data_obj$S)
  }
  data_obj$min_det_matrix <- diag(1, data_obj$n_det + 1, data_obj$n_det) # matrix used for .p 
  
  ## Create start parameter values
  # density_pars <- rep(0, ncol(model.matrix(data_obj$gam_fit)))
  density_pars <- c(-12, 45, -53)
  names(density_pars) <- colnames(model.matrix(data_obj$gam_fit))
  
  parameters <- c(logit_g0 = 0 , 
                  log_beta_r = log(19) , 
                  log_sd_r = log(2.7) , 
                  log_mu_s = log(163))
  if (!FIXED_SL) parameters <- c(parameters,
                                 log_sd_s = log(5))
  if (USE_BEARINGS == 1) parameters <- c(parameters,
                                         log_kappa = log(5))
  if (USE_BEARINGS == 2) parameters <- c(parameters,
                                         log_kappa_low = log(0.3) ,
                                         log_kappa_high = log(37) ,
                                         logit_mix_bear = -0.5)
  parameters <- c(parameters, density_pars)
  
  data_obj$dens_par_names <- colnames(model.matrix(data_obj$gam_fit))
  
  ## ===========================================================================
  ## 2.2 SOURCE OTHER R AND CPP SCRIPTS
  ## ---------------------------------------------------------------------------
  source("Scripts/bowhead whales/nllR.R")
  ## ---------------------------------------------------------------------------
  
  ## ===========================================================================
  ## 2.3 RUN THE OPTIMISER 
  ## ---------------------------------------------------------------------------
  res <- nlminb(start = parameters, objective = nllR, dat = data_obj,
                # lower = -50,
                # upper = 50,
                control = list(trace = 1))
  ## ---------------------------------------------------------------------------
  
  ## ===========================================================================
  ## 2.4 ADD OTHER INFORMATION TO DATA OBJECT, AND RETURN
  ## ---------------------------------------------------------------------------
  fit_obj <- list(dat = data_obj,
                  start_pars = parameters,
                  fit_res = res)
  real_par <- exp(res$par)
  real_par[which(names(res$par) %in% c("logit_g0", "logit_mix_par"))] <- 
    exp(res$par[which(names(res$par) %in% c("logit_g0", "logit_mix_par"))]) /
    (1 + exp(res$par[which(names(res$par) %in% c("logit_g0", "logit_mix_par"))]))
  fit_obj$real_par <- real_par
  fit_obj$aic <- -2 * res$value + 2 * length(res$par) # Burnham & Anderson (2002)
  fit_obj$aicc <- fit_obj$aic + 2 * length(res$par) * (length(res$par) + 1) / 
    (data_obj$n_call - length(res$par) - 1) # Burnham & Anderson (2002)
  fit_obj$bic <- -2 * res$value + length(res$par) * log(data_obj$n_call)
  
  gam_fit <- data_obj$gam_fit
  gam_fit$coefficients <- res$par[names(gam_fit$coefficients)]
  fit_obj$D <- data.frame(area = data_obj$A_x, 
                          density = exp(predict(gam_fit, 
                                                data = model.matrix(gam_fit))))
  fit_obj$D$N <- fit_obj$D$area * fit_obj$D$density
  fit_obj$N_total <- sum(fit_obj$D$N)
  
  ## Return the fit object
  return(fit_obj)
  ## ---------------------------------------------------------------------------
}) 

## -----------------------------------------------------------------------------








## =============================================================================
## 
## -----------------------------------------------------------------------------



## -----------------------------------------------------------------------------