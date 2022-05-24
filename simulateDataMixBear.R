# rm(list = setdiff(ls(), c("f_density", "par", "cov_density", "detectors",
#                           "min_no_detections", "det_function",
#                           "SINGLE_SL", "WITH_NOISE")))

simulateDataMixBear <- function(par, f_density, cov_density, detectors,
                         min_no_detections, det_function, SINGLE_SL = FALSE, 
                         WITH_NOISE = FALSE, MIX_BEAR = FALSE, ...) {
  # Description:
  #   Based on the inputs, creates a density of the area. Based on this, 
  #   capture histories are simulated and returned. This can be used to test the
  #   accuracy of a spatial capture-recapture fitting method.
  
  # Inputs:
  #   
  
  # Outputs:
  #   dat         - [list] 
  
  
  #----------------# Load libraries and perform input checks #-----------------#
  library(dplyr)
  library(mgcv)
  library(matrixStats)
  library(circular)
  library(truncnorm)
  library(raster)
  library(geosphere)
  
  source("Scripts/Bowhead Whales/hidden_functions.R")
  source("Scripts/Bowhead Whales/plotDensity.R")
  
  #-------------------------# Input checks #-----------------------------------#
  
  #-------------------------# Extract data #-----------------------------------#
  n_det <- nrow(detectors) # number of detectors
  n_grid <- nrow(cov_density) # number of grid points
  
  args <- list(...)

  trunc_level <- args$trunc_level
  # trunc_level <- 75
  
  # Extract parameters and convert to the real scale (add errors to ensure 
  # correct domain)
  par_det <- par$par_det
  
  if (MIX_BEAR) {
    kappa_low <- par$par_bear["kappa_low"]
    kappa_high <- par$par_bear["kappa_high"]
    mix_par <- par$par_bear["mix_par"]
  } else {
    kappa <- par$par_bear["kappa"] 
  }
  par_rl <- par$par_rl
  par_sl <- par$par_sl
  par_dens <- par$par_dens # On the log scale, as D = exp(...)
  if (WITH_NOISE) {par_noise <- par$par_noise}
  
  if (WITH_NOISE) {
    if (SINGLE_SL) {
      cat(paste0("Simulating data with:\n",
                 "  A single source level of ", par$par_sl["mu"], ",\n",
                 "  and a signal-to-noise ratio trunctation value of ", trunc_level, ".\n"))
    } else {
      cat(paste0("Simulating data with:\n",
                 "  A variables source level ~ N(", par$par_sl["mu"], ", ", 
                 par$par_sl["sd"], "),\n",
                 "  and a signal-to-noise ratio trunctation value of ", trunc_level, ".\n"))
    }
  } else {
    if (SINGLE_SL) {
      cat(paste0("Simulating data with:\n",
                 "  A single source level of ", par$par_sl["mu"], ",\n",
                 "  and a received level trunctation value of ", trunc_level, ".\n"))
    } else {
      cat(paste0("Simulating data with:\n",
                 "  A variables source level ~ N(", par$par_sl["mu"], ", ", 
                 par$par_sl["sd"], "),\n",
                 "  and a received level trunctation value of ", trunc_level, ".\n"))
    }
  }
  
  # Extract density function parameters
  vars <- all.vars(f_density)
  vars <- c("(Intercept)", vars[-1])
  
  # Check if all density parameters have starting values
  if (!all(vars %in% names(par_dens))) {
    stop("Not all density parameters have starting values.")
  }
  
  # Extract the area
  A <- subset(cov_density, select = c(long, lat, area))
  
  # Potentially scale the grid_density data
  cov_density_scaled <- subset(cov_density, select = -area)
  cov_density_scaled[, -c(1, 2)] <- scale(subset(cov_density_scaled, 
                                                 select = -c(long, lat)))
  
  # Create a matrix of distances
  distances <- pointDistance(p1 = cov_density_scaled[, c("long", "lat")], 
                             p2 = detectors[, c("long", "lat")], 
                             lonlat = TRUE)
  
  #----------------------# Start simulating data #-----------------------------#
  
  # Create a gam object
  gam_fit <- gam(f_density, data = cbind(D = 0, cov_density_scaled))
  if (length(gam_fit$smooth) != 0) {
    k <- gam_fit$smooth[[1]]$bs.dim
    gam_par <- rep(0, k - 1)
    smooth_terms <- paste0(gam_fit$smooth[[1]]$label, ".", 1:(k - 1))
    names(gam_par) <- smooth_terms
    
    par <- c(par, gam_par)
  } else {smooth_terms <- NULL}
  
  # Extract design matrix
  design <- model.matrix(gam_fit)
  
  # Derive the density
  D <- .densityGAM(x = design, gam_fit = gam_fit, par = par_dens)
  if (any(D < 0)) {
    print("At least one element in D is negative!")
  }
  # Print map of the density
  map <- plotDensity(d = data.frame(long = cov_density$long, 
                                    lat = cov_density$lat, 
                                    density = D))
  print(map)
  
  # set.seed(1)
  
  ## Determine how many calls were simulated and extract distances to detectors
  emitted_call_per_cell <- rpois(n = n_grid, lambda = A$area * D)
  n_call <- sum(emitted_call_per_cell)
  cat(paste0(n_call, " bowhead whale calls were simulated.\n"))
  
  call_sample_index <- rep(1:n_grid, times = emitted_call_per_cell)
  distances <- distances[call_sample_index, ]
  calls <- cov_density[call_sample_index, ]
  
  # set.seed(1)
  ## Simulate sound characteristics for every call
  if (SINGLE_SL) {
    source_levels <- rep(par_sl["mu"], n_call) 
  } else {
    source_levels <- rnorm(n_call, par_sl["mu"], par_sl["sd"])
    # source_levels <- rtruncnorm(n_call,
    #                             a = par_sl["lower"],
    #                             b = par_sl["upper"],
    #                             mean = par_sl["mu"],
    #                             sd = par_sl["sd"])
  }
  
  # For every call, derive received level and add measurement error
  received_levels <- matrix(source_levels, nrow = n_call, ncol = n_det,
                            byrow = FALSE) -
    par_rl["beta"] * log(distances, base = 10) + # transmission loss
    matrix(rnorm(n_call * n_det, 0, par_rl["sigma"]), nrow = n_call) # measurement error

  if (WITH_NOISE) {
    # Simulate mean noise level for every call
    noise <- rtruncnorm(n_call,
                        a = par_noise["lower"],
                        b = par_noise["upper"],
                        mean = par_noise["mu"],
                        sd = par_noise["sd"])
    # For every noise level, create slight deviations for all detectors
    noise <- matrix(rep(noise, each = n_det), nrow = n_call, ncol = n_det,
                    byrow = TRUE)
    # noise <- noise + matrix(rnorm(length(noise), 0, 2), nrow = n_call, #hard-coded sd = 2
    #                         ncol = n_det, byrow = TRUE)
    
    # Derive the signal to noise ratio
    snr <- received_levels - noise
  }

  
  #----------------# Create detection histories #------------------------------#
  sd_r <- par_rl["sigma"]

  # Derive detection probabilities based on snr 
  # if (det_function == "half-normal") {
  #   det_probs <- .gHN(distances = distances, par = par_det)
  # } else if (WITH_NOISE) {
  #   det_probs <- .gSNR(snr = snr, par = par_det, type = det_function, 
  #                      sd_r = sd_r, trunc_level = trunc_level)
  # } else {
  #   det_probs <- .gSNR(snr = rl, par = par_det, type = det_function, 
  #                      sd_r = sd_r, trunc_level = trunc_level)
  # }
  
  if (det_function == "simple") {
    if (WITH_NOISE) {
      det_probs <- par_det["g0"] * (snr > trunc_level)
    } else {
    det_probs <- par_det["g0"] * (received_levels > trunc_level)
    } 
  }
    
  # OPTION 2 : Add rl error AFTER det probs are derived
  # received_levels <- received_levels +
  #   matrix(rnorm(n_call * n_det, 0, par_rl["sigma"]), nrow = n_call) # measurement error
  
  # set.seed(1)
  # Create detection history consisting of 1's and 0's
  det_hist <- det_probs > runif(n_call * n_det, min = 0, max = 1)
  mode(det_hist) <- "numeric"
  
  # Create index for whether min_no_detections was satisfied
  ENOUGH_DETS <- rowSums(det_hist) >= min_no_detections
  detected_calls <- calls[ENOUGH_DETS, ]
  
  # Print result of simulation
  cat(paste0(sum(ENOUGH_DETS), " bowhead whale calls were detected at least ", 
             min_no_detections, " times.\n"))    
  
  #---------------------# Create bearing histories #---------------------------#
  # Derive bearings for all calls with enough detections
  bearings <- t(apply(detected_calls[, c("long", "lat")], 1, function(grid_point) {
    bear <- geosphere::bearing(p1 = as.matrix(detectors[, c("long", "lat")]),
                               p2 = grid_point)
  }))
  bearings <- circular(x = bearings, template = "geographics", modulo = "2pi",
                       units = "degrees")
  # set.seed(1)
  # Add the measurement error
  if (MIX_BEAR) {
    tests <- matrix(runif(length(bearings), 0, 1), 
                    ncol = n_det, 
                    nrow = nrow(bearings))
    tests <- tests < mix_par
    
    high_errors <- matrix(circular::rvonmises(n = length(bearings), 
                                              mu = circular(0, template = "geographics", 
                                                            modulo = "2pi", units = "degrees"), 
                                              kappa = kappa_low), 
                          ncol = n_det, nrow = nrow(bearings))
    
    low_errors <- matrix(circular::rvonmises(n = length(bearings), 
                                             mu = circular(0, template = "geographics", 
                                                           modulo = "2pi", units = "degrees"), 
                                             kappa = kappa_low + kappa_high), 
                         ncol = n_det, nrow = nrow(bearings))
    
    errors <- low_errors
    errors[tests] <- high_errors[tests]
    # errors <- mix_par * 
    #   matrix(circular::rvonmises(n = length(bearings), 
    #                              mu = circular(0, template = "geographics", 
    #                                            modulo = "2pi", units = "degrees"), 
    #                            kappa = kappa_low), 
    #          ncol = n_det, nrow = nrow(bearings)) + (1 - mix_par) *
    #   matrix(circular::rvonmises(n = length(bearings), 
    #                              mu = circular(0, template = "geographics", 
    #                                            modulo = "2pi", units = "degrees"), 
    #                              kappa = kappa_low + kappa_high), 
    #          ncol = n_det, nrow = nrow(bearings))
  } else {
    errors <- matrix(circular::rvonmises(n = length(bearings), 
                                         mu = circular(0, template = "geographics", 
                                                       modulo = "2pi", units = "degrees"), 
                                         kappa = kappa),
                     ncol = n_det, nrow = nrow(bearings))
  }
  bearings <- bearings + errors
  bearings <- circular(x = bearings, template = "geographics", 
                       modulo = "2pi", units = "degrees")

  ################# Add the data for detected calls to dat #####################
  
  ## Create output data
  det_hist <- det_hist[ENOUGH_DETS, ]
  NO_DETECTION <- !det_hist
  
  received_level_hist <- received_levels[ENOUGH_DETS, ]
  received_level_hist[NO_DETECTION] <- NA

  if (WITH_NOISE) {
    noise_call <- noise[ENOUGH_DETS, ]
    noise_random <- noise[sample(1:n_call, 1000, replace = TRUE), ]
    noise_random <- rtruncnorm(1000,
                               a = par_noise["lower"],
                               b = par_noise["upper"],
                               mean = par_noise["mu"],
                               sd = par_noise["sd"])
    noise_random <- matrix(rep(noise_random, each = n_det), ncol = n_det, byrow = TRUE)
    
    # noise_random <- noise_random + matrix(rnorm(length(noise_random), 0, 2), nrow = 1000,  # hardcoded sd = 2
    #                                       ncol = n_det, byrow = TRUE)
  }
  
  bearing_hist <- bearings
  bearing_hist[NO_DETECTION] <- NA
  
  # Only essential added right now, could output more info but not required.
  if (WITH_NOISE) {
    dat <- list(det_hist = det_hist,
                received_levels = received_level_hist,
                noise_call = noise_call,
                noise_random = noise_random,
                bearings = bearing_hist)
  } else {
    dat <- list(det_hist = det_hist,
                received_levels = received_level_hist,
                bearings = bearing_hist)
  }

  return(dat)
}
  