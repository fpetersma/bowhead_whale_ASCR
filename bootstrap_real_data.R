## =============================================================================
## Bootstrap the real data to 999 data sets
## 
## Felix Petersma
## =============================================================================

## =============================================================================
## 1. LOAD DATA AND LIBRARIES
## -----------------------------------------------------------------------------

## Libraries
library(circular)
library(tidyverse)
library(Rcpp)       # for combining R and cpp
library(parallel)   # for parallel processing
library(mgcv)

## Load best model
load("Real data output/31_fits_real_data_all_convergence.RData")
best_model <- fits[[31]]

dat <- best_model$dat
dat$S

set.seed(18092021)

## Resample with replacement
bootstrapped_data <- lapply(1:999, function(i) {
  new_indices <- sample(x = 1:dat$n_call, 
                        size = dat$n_call, 
                        replace = TRUE)
  boot <- dat
  boot$Y_rec <- boot$Y_rec[new_indices, ]
  boot$W <- boot$W[new_indices, ]
  boot$R <- boot$R[new_indices, ]
  
  return(boot)
})

## Create start values
# density_pars <- best_model$fit_res$par[-(1:8)]
parameters <- best_model$fit_res$par

## =============================================================================
## 2. SOURCE OTHER R AND CPP SCRIPTS, AND FIT TO BOOTSTRAPPED DATA
## -----------------------------------------------------------------------------
source("Scripts/bowhead whales/nllR.R")

results <- list()

for (i in 1:100) {
  ## Fit model
  res <- nlminb(start = parameters, 
                objective = nllR, 
                dat = bootstrapped_data[[i]],
                # lower = -50,
                # upper = 50,
                control = list(trace = 1, 
                               eval.max = 1000, 
                               iter.max = 500))
  
  ## Add other info to fit object
  fit_obj <- list(dat = data_obj,
                  start_pars = parameters,
                  fit_res = res)
  real_par <- exp(res$par)
  real_par[which(names(res$par) %in% c("logit_g0", "logit_mix_par"))] <- 
    exp(res$par[which(names(res$par) %in% c("logit_g0", "logit_mix_par"))]) /
    (1 + exp(res$par[which(names(res$par) %in% c("logit_g0", "logit_mix_par"))]))
  fit_obj$real_par <- real_par
  ## res$objective is the negative log likelihood, thus the minusses cancel out
  fit_obj$aic <- 2 * res$objective + 2 * length(res$par) # Burnham & Anderson (2002)
  fit_obj$aicc <- fit_obj$aic + 2 * length(res$par) * (length(res$par) + 1) / 
    (data_obj$n_call - length(res$par) - 1) # Burnham & Anderson (2002)
  fit_obj$bic <- 2 * res$objective + length(res$par) * log(data_obj$n_call)
  
  gam_fit <- data_obj$gam_fit
  gam_fit$coefficients <- res$par[names(gam_fit$coefficients)]
  fit_obj$D <- data.frame(area = data_obj$A_x, 
                          density = exp(predict(gam_fit, 
                                                data = model.matrix(gam_fit))))
  fit_obj$D$N <- fit_obj$D$area * fit_obj$D$density
  fit_obj$N_total <- sum(fit_obj$D$N)
  
  results[[i]] <- fit_obj
}

## =============================================================================
## CREATE FINAL TABLE
## -----------------------------------------------------------------------------

fit_overview <- data.frame(ID = 1:length(fits), 
                           density_model = sapply(fits, function(fit) Reduce(paste, deparse(fit$dat$f_density))),
                           N_hat = sapply(fits, function(fit) fit$N_total),
                           convergence = sapply(fits, function(fit) fit$fit_res$convergence),
                           n_par = sapply(fits, function(fit) length(fit$fit_res$par)),
                           llk = sapply(fits, function(fit) fit$fit_res$objective),
                           AIC = sapply(fits, function(fit) fit$aic),
                           AICc = sapply(fits, function(fit) fit$aicc),
                           BIC = sapply(fits, function(fit) fit$bic))

fit_overview <- fit_overview[order(fit_overview$AIC, decreasing = TRUE), ]
## =============================================================================
