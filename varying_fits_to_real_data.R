## =============================================================================
## Fitting models to real data
##
## Felix Petersma
##
## Created on 02/09/2021
##
## Update 18/11/2022: potentially have to do a rerun with on localised calls.
## =============================================================================


###### Create models and start values ################################
models <- list(D ~ 1,
               D ~ distance_to_coast + distance_to_coast2,
               D ~ distance_to_coast + distance_to_coast2 + distance_to_coast3,
               D ~ depth,
               D ~ depth + depth2,
               D ~ logdepth,
               D ~ logdepth + depth + depth2,
               D ~ logdepth + distance_to_coast + distance_to_coast2 + distance_to_coast3,
               D ~ logdepth + depth + distance_to_coast + distance_to_coast2,
               D ~ logdepth + depth + depth2 + distance_to_coast,
               D ~ depth + distance_to_coast,
               D ~ depth + distance_to_coast + distance_to_coast2,
               D ~ depth + depth2 + distance_to_coast,
               D ~ depth + depth2 + distance_to_coast + distance_to_coast2,
               D ~ depth + depth2 + distance_to_coast + distance_to_coast2 + distance_to_coast3,
               D ~ s(depth, k = 3, fx = TRUE),
               D ~ s(depth, k = 4, fx = TRUE),
               D ~ s(depth, k = 5, fx = TRUE),
               D ~ s(depth, k = 6, fx = TRUE),
               D ~ s(depth, k = 7, fx = TRUE),
               D ~ s(depth, k = 8, fx = TRUE),
               D ~ s(depth, k = 6, fx = TRUE) + distance_to_coast,
               D ~ s(depth, k = 6, fx = TRUE) + distance_to_coast + distance_to_coast2,
               D ~ s(distance_to_coast, k = 3, fx = TRUE),
               D ~ s(distance_to_coast, k = 4, fx = TRUE),
               D ~ s(distance_to_coast, k = 5, fx = TRUE),
               D ~ s(distance_to_coast, k = 6, fx = TRUE),
               D ~ s(distance_to_coast, k = 7, fx = TRUE),
               D ~ s(distance_to_coast, k = 8, fx = TRUE),
               D ~ s(distance_to_coast, k = 6, fx = TRUE) + depth,
               D ~ s(distance_to_coast, k = 6, fx = TRUE) + depth + depth2,
               D ~ s(depth, k = 4, fx = TRUE) + s(distance_to_coast, k = 4, fx = TRUE),
               D ~ s(depth, k = 6, fx = TRUE) + s(distance_to_coast, k = 6, fx = TRUE),
               D ~ distance_to_coast + depth + depth:distance_to_coast,
               D ~ distance_to_coast + distance_to_coast2 + depth + depth2 + depth:distance_to_coast)

# models <- list(D ~ distance_to_coast + depth + depth:distance_to_coast)



## =============================================================================
## 1. LOAD DATA AND LIBRARIES
## -----------------------------------------------------------------------------

## Libraries
library(circular)
library(tidyverse)
library(Rcpp)       # for combining R and cpp
library(parallel)   # for parallel processing
library(mgcv)       # for smooths
library(pbapply)

## Constants
A_s <- 3
FIXED_SL <- TRUE
USE_BEARINGS <- 2
WITH_NOISE <- FALSE
trunc_level <- 96


## Load the real data
detections <- read.csv("Data/detections_31-08-2010_successful-loc.csv")
bearings <- read.csv("Data/bearings_31-08-2010_successful-loc.csv")
received_levels <- read.csv("Data/received_levels_31-08-2010_successful-loc.csv")

DASAR <- as.data.frame(readr::read_tsv("Data/DASARs.txt"))

mesh <- read.csv("Data/alaska_albers_grid_adaptive_levels=2_inner=10k_outer=50k_maxD2C=Inf_area=8450_n=438.csv")
# 
# fits <- list()
# for (model_index in 22) {

## Create parallel processing
cl <- makeCluster(21)
clusterExport(cl, varlist = ls())

model_index <- 8

fits <- pblapply(1:21, function(model_index) {
  # model_index <- 33
  cat("Evaluating model:", model_index, "\n")
  ## Libraries
  library(circular)   # for von Mises distribution
  library(mgcv)       # for smooths
  library(dplyr)      # for piping
  
  ## Extract the current density  model specification
  f_density <- models[[model_index]]
  
  ## Create data object to store data in
  dat <- list()
  
  ## Subset DASAR detector data
  DASAR$long <- -abs(DASAR$long) # Turn longitude from westing to easting (more common)
  detectors <- DASAR[DASAR$year == 2010 & DASAR$site == 5, c("long", "lat")] 
  
  ## Create standardised covariates data.frame
  covariates <- dplyr::select(mesh, c(distance_to_coast, depth)) %>% 
    mutate(distance_to_coast2 = distance_to_coast ^ 2, 
           distance_to_coast3 = distance_to_coast ^ 3,
           depth = abs(depth),
           depth2 = depth ^ 2,
           logdepth = log(depth)) 
  dat$covariates <- scale(covariates, center = FALSE, scale = TRUE)
  
  dat$A_x <- mesh$area
  
  ## Run next section to filter on received levels
  if (WITH_NOISE) {
    sufficient_snr <- received_levels - noise_call >= trunc_level # get indices 
    sufficient_snr[is.na(sufficient_snr)] <- FALSE
    
    detections <- sufficient_snr * 1
    enough_dets <- Rfast::rowsums(detections) > 1
    
    # only keep calls with enough detections, after truncation at 15dB snr
    sufficient_snr <- sufficient_snr[enough_dets, ]
    detections <- detections[enough_dets, ]
    
    bearings <- as.matrix(bearings)[enough_dets, ]
    bearings[!sufficient_snr] <- NA
    
    noise_call <- as.matrix(noise_call)[enough_dets, ]
    
    received_levels <- as.matrix(received_levels)[enough_dets, ]
    received_levels[!sufficient_snr] <- NA
  } else if (!WITH_NOISE) {
    sufficient_rl <- received_levels  >= trunc_level # get indices 
    sufficient_rl[is.na(sufficient_rl)] <- FALSE
    
    detections <- sufficient_rl * 1
    enough_dets <- Rfast::rowsums(detections) > 1
    
    # only keep calls with enough detections, after truncation at 15dB snr
    sufficient_rl <- sufficient_rl[enough_dets, ]
    detections <- detections[enough_dets, ]
    
    bearings <- as.matrix(bearings)[enough_dets, ]
    bearings[!sufficient_rl] <- NA
    
    # noise_call <- as.matrix(noise_call)[enough_dets, ]
    
    received_levels <- as.matrix(received_levels)[enough_dets, ]
    received_levels[!sufficient_rl] <- NA
  }
  
  # index <- sample(1:nrow(detections), sample_size)
  det_hist <- detections#[index, ]
  bearings_hist <- bearings#[index, ]
  received_levels_hist <- received_levels#[index, ]
  
  dat$det_hist <- det_hist
  dat$received_levels <- received_levels
  
  # Create a matrix of distances
  distances <- raster::pointDistance(p1 = mesh[, c("long", "lat")], 
                             p2 = detectors[, c("long", "lat")], 
                             lonlat = TRUE)
  
  
  bearings_deg <- circular(bearings,
                           units = "degrees",
                           template = "geographics") # 'geographics' means clockwise rotation with 0 at north
  bearings_rad <- conversion.circular(bearings_deg)
  dat[["bearings_rad"]] <- bearings_rad # Add bearings as circular in radians to dat
  
  # Get degrees bearings from detectors to all grid points 
  grid_bearings <- t(apply(mesh[, c("long", "lat")], 1, function(grid_point) {
    bear <- geosphere::bearing(p1 = as.matrix(detectors[, c("long", "lat")]),
                               p2 = grid_point)
  }))
  # First convert to 'circular' 
  grid_bearings <- circular(grid_bearings, template = "geographics", 
                            modulo = "2pi", units = "degrees")
  # Then convert to radians
  grid_bearings <- conversion.circular(grid_bearings, units = "radians")
  dat[["grid_bearings"]] <- grid_bearings
  
  dat$distances <- distances
  dat$bearings_rad <- bearings_rad
  
  dat$trunc_level <- trunc_level
  
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
                   # A_s = A_s,
                   trunc_level = dat$trunc_level,
                   # design_matrix = design_matrix,
                   gam_fit = gam_fit,
                   # S = seq(from = 80 + A_s / 2, to = 200 - A_s / 2, by = A_s),
                   n_det = ncol(dat$distances),
                   n_grid = nrow(dat$distances),
                   n_call = nrow(dat$det_hist),
                   FIXED_SL = FIXED_SL,
                   USE_BEARINGS = USE_BEARINGS) # create data
  if (!FIXED_SL) {
    data_obj$A_s <- A_s
    data_obj$S <- seq(from = 100 + A_s / 2, to = 220 - A_s / 2, by = A_s)
    data_obj$n_sl <- length(data_obj$S)
  }
  data_obj$min_det_matrix <- diag(1, data_obj$n_det + 1, data_obj$n_det) # matrix used for .p 
  
  ## Create start parameter value
  density_pars <- rep(0, ncol(model.matrix(data_obj$gam_fit)))
  names(density_pars) <- colnames(model.matrix(data_obj$gam_fit))
  
  parameters <- c(logit_g0 = 0 , 
                  log_beta_r = log(18) , 
                  log_sd_r = log(2) , 
                  log_mu_s = log(163))
  if (!FIXED_SL) parameters <- c(parameters,
                                 log_sd_s = 2)
  if (USE_BEARINGS == 1) parameters <- c(parameters,
                                         log_kappa = 1)
  if (USE_BEARINGS == 2) parameters <- c(parameters,
                                         log_kappa_low = 1 ,
                                         log_kappa_high = log(35) ,
                                         logit_mix_bear = -0.5)
  parameters <- c(parameters, density_pars)
  
  data_obj$dens_par_names <- colnames(model.matrix(data_obj$gam_fit))
  ## =============================================================================
  
  ## =============================================================================
  ## 2. SOURCE OTHER R AND CPP SCRIPTS
  ## -----------------------------------------------------------------------------
  source("Scripts/bowhead whales/nllR.R")
  
  ## =============================================================================

  ## =============================================================================
  ## 3. RUN THE OPTIMISER 
  ## -----------------------------------------------------------------------------
  # system.time({
  # res <- optim(par = parameters, fn = nllR, dat = data_obj, method = "Nelder-Mead",
  #              hessian = FALSE,
  #              # lower = -50,
  #              # upper = 50,
  #              control = list(trace = 6, 
  #                             REPORT = TRUE,
  #                             maxit = 2000))
  res <- nlminb(start = parameters, objective = nllR, dat = data_obj,
                # lower = -50,
                # upper = 50,
                control = list(trace = 1,
                               rel.tol = 1e-8,
                               eval.max = 1000,
                               iter.max = 500))
  # })
  ## =============================================================================
  #  Truth:             0.405     2.701   0.693     5.011   0.000       2.944   -2.197   -1        2       -1.5
  #  35:     4263.9473: 0.452523  2.71824 0.700802  5.01431 0.00945608  2.82123 -2.38982 -1.13384  2.48271 -2.22861
  #  33:     4300.1762: 0.453393  2.72198 0.718492  5.01494 -0.198829   2.93911 -2.20946 -1.20934  2.86865 -2.69616
  #  34:     4226.5505: 0.477780  2.72429 0.673404  5.01761 -0.177191   2.90168 -2.35719 -1.07825  2.10859 -1.71153
  ## =============================================================================
  ## 4. ADD OTHER INFORMATION TO DATA OBJECT
  ## -----------------------------------------------------------------------------
  fit_obj <- list(dat = data_obj,
                  start_pars = parameters,
                  fit_res = res)
  real_par <- exp(res$par)
  real_par[which(names(res$par) %in% c("logit_g0", "logit_mix_par"))] <- 
    exp(res$par[which(names(res$par) %in% c("logit_g0", "logit_mix_par"))]) /
    (1 + exp(res$par[which(names(res$par) %in% c("logit_g0", "logit_mix_par"))]))
  fit_obj$real_par <- real_par
  fit_obj$aic <- 2 * res$objective + 2 * length(res$par) # Burnham & Anderson (2002)
  fit_obj$aicc <- fit_obj$aic + 2 * length(res$par) * (length(res$par) + 1) / 
    (data_obj$n_call - length(res$par) - 1) # Burnham & Anderson (2002)
  fit_obj$bic <- 2 * res$objective + length(res$par) * log(data_obj$n_call)
  
  gam_fit <- data_obj$gam_fit
  gam_fit$coefficients <- res$par[names(gam_fit$coefficients)]
  fit_obj$D <- data.frame(area = data_obj$A_x, 
                          density = exp(predict(gam_fit, data = model.matrix(gam_fit))))
  fit_obj$D$N <- fit_obj$D$area * fit_obj$D$density
  fit_obj$N_total <- sum(fit_obj$D$N)
  fit_obj$par_hist <- read.csv(file = paste0("parameter_history_Rcpp_", 
                                             f_density[3], ".csv"))

  ## =============================================================================
  ## Below no longer need for pbapply version
  # fits[[model_index]] <- fit_obj
  return(fit_obj)
}, cl = cl); parallel::stopCluster(cl);

## =============================================================================
## CREATE FINAL TABLE
## -----------------------------------------------------------------------------

fit_overview <- data.frame(ID = 1:length(fits), 
                           density_model = sapply(fits, function(fit) Reduce(paste, deparse(fit$dat$f_density))),
                           N_hat = sapply(fits, function(fit) fit$N_total),
                           convergence = sapply(fits, function(fit) fit$fit_res$convergence),
                           n_par = sapply(fits, function(fit) length(fit$fit_res$par)),
                           llk = sapply(fits, function(fit) -fit$fit_res$objective),
                           AIC = sapply(fits, function(fit) fit$aic),
                           AICc = sapply(fits, function(fit) fit$aicc),
                           BIC = sapply(fits, function(fit) fit$bic))

fit_overview <- fit_overview[order(fit_overview$AIC), ]
## =============================================================================

## Extract hessians if required [added on 01/12/22 after suggestion by Len]
source("Scripts/bowhead whales/nllR.R")
hessian <- numDeriv::hessian(func = nllR,
                             x = fits[[8]]$fit_res$par,
                             dat = fits[[8]]$dat)
vcov <- solve(-hessian, tol = 3e-49)

# alternative to numDeriv... pracma!
hessian <- pracma::hessian(f = nllR,
                           x0 = fits[[33]]$fit_res$par,
                           dat = fits[[33]]$dat)

vcov <- solve(-hessian)

