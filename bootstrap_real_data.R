## =============================================================================
## Bootstrap the real data to 999 data sets
## 
## Felix Petersma
##
## [09/12/2022] updated for model 33.
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
library(pbapply)

## Load best model
# load("fits_1_21_nlminb_reltol=1e-8_n=443.RData") # Load list of fits
load("fits_1_35_nlminb_n=443.RData")
best_model <- fits[[33]]

dat <- best_model$dat

set.seed(17051993)

## Resample the calls with replacement
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
start_pars <- best_model$fit_res$par
# start_pars[9:13] <- start_pars[9:13] / 100 # divide by 100 for improved convergence
# start_pars["(Intercept)"] <- -11

## =============================================================================
## 2. SOURCE OTHER R AND CPP SCRIPTS, AND FIT TO BOOTSTRAPPED DATA
## -----------------------------------------------------------------------------
source("Scripts/Bowhead Whales/nllR.R") # no parallel here.

# nllR(start_pars, bootstrapped_data[[1]])
rm(best_model, dat, fits)
# results <- list()
start <- 1
end <- 999

cl <- makeCluster(30, outfile="bootstrapping_console_output.txt")
clusterExport(cl, c("start_pars", "bootstrapped_data", "nllR"))

# for (i in start:end) {
results <- pblapply(start:end, function(j) { # parallel here!
  cat("Currently at run: ", j, "\n", sep = "")
  dat_j <- bootstrapped_data[[j]]
  ## Fit model
  res <- nlminb(start = start_pars, 
                objective = nllR, 
                dat = dat_j,
                # lower = -50,
                # upper = 50,
                control = list(trace = 0, 
                               # rel.tol = 1e-8,
                               eval.max = 1000, 
                               iter.max = 500))
  
  # save(list = "res", file = paste0("test_results_", j, ".RData"))
  
  ## Add other info to fit object
  fit_obj <- list(dat = dat_j,
                  start_pars = start_pars,
                  fit_res = res)
  real_par <- exp(res$par)
  real_par[which(names(res$par) %in% c("logit_g0", "logit_mix_par"))] <-
    exp(res$par[which(names(res$par) %in% c("logit_g0", "logit_mix_par"))]) /
    (1 + exp(res$par[which(names(res$par) %in% c("logit_g0", "logit_mix_par"))]))
  fit_obj$real_par <- real_par
  ## res$objective is the negative log likelihood, thus the minusses cancel out
  fit_obj$aic <- 2 * res$objective + 2 * length(res$par) # Burnham & Anderson (2002)
  fit_obj$aicc <- fit_obj$aic + 2 * length(res$par) * (length(res$par) + 1) /
    (dat_j$n_call - length(res$par) - 1) # Burnham & Anderson (2002)
  fit_obj$bic <- 2 * res$objective + length(res$par) * log(dat_j$n_call)

  gam_fit <-  dat_j$gam_fit
  gam_fit$coefficients <- res$par[names(gam_fit$coefficients)]
  fit_obj$D <- data.frame(area =  dat_j$A_x,
                          density = exp(mgcv::predict.gam(gam_fit,
                                                data = model.matrix(gam_fit))))
  fit_obj$D$N <- fit_obj$D$area * fit_obj$D$density
  fit_obj$N_total <- sum(fit_obj$D$N)

  return(fit_obj)
  
  # results[[i]] <- fit_obj
  # save(list = "results", file = paste0("results_", start, "-", end, ".RData"))
}, cl = cl); stopCluster(cl)
save(list = "results", file = "results_model8_750_999_.RData")
# save.image("first_200_bootstraps.RData")