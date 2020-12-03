# rm(list = setdiff(ls(), c("f_density", "par", "cov_density", "detectors",
#                           "A", "min_no_detections", "det_function")))

simulateData <- function(par, f_density, cov_density, detectors, A,
                         min_no_detections, det_function) {
  # Description:
  #   Based on the inputs, creates a density of the area. Based on this, 
  #   capture histories are simulated and returned. This can be used to test the
  #   accuracy of a spatial capture-recapture fitting method.
  
  # Inputs:
  #   
  
  # Outputs:
  #   dat         - [list] 
  
  
  ################## Load libraries and perform input checks ###################
  library(dplyr)
  library(mgcv)
  library(matrixStats)
  library(circular)
  library(truncnorm)
  
  source("Scripts/Bowhead Whales/hidden_functions.R")
  source("Scripts/Bowhead Whales/plotDensity.R")
  
  ########################### Input checks #####################################
  
  ########################### Extract data #####################################
  n_det <- nrow(detectors) # number of detectors
  n_grid <- nrow(cov_density) # number of grid points
  
  # Extract parameters and convert to the real scale (add errors to ensure 
  # correct domain)
  par_det <- par$par_det
  kappa <- par$par_bear["kappa"] 
  par_rl <- par$par_rl
  par_sl <- par$par_sl
  par_noise <- par$par_noise
  par_dens <- par$par_dens # On the log scale, as D = exp(...)
  
  # Extract density function parameters
  vars <- all.vars(f_density)
  vars <- c("(Intercept)", vars[-1])
  
  # Check if all density parameters have starting values
  if (!all(vars %in% names(par_dens))) {
    stop("Not all density parameters have starting values.")
  }
  
  ## Potentially scale the grid_density data
  cov_density_scaled <- cov_density
  cov_density_scaled[, -c(1, 2)] <- scale(select(cov_density, -c("x", "y")))
  
  ## Create a matrix of distances
  distances <- apply(detectors, 1, function(z) {
    y <- t(cov_density_scaled[, c("x", "y")]) - z
    y <- y ^ 2
    return(sqrt(colSums(y)))
  })
  
  ######################## Start simulating data ###############################
  
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
  map <- plotDensity(d = data.frame(x = cov_density$x, y = cov_density$y, density = D))
  print(map)
  
  # create detection function plot
  x <- 0:40
  plot(x = x, 
       y = .gSNR(snr = x, par = par_det, type = det_function),
       main = "Detection function", 
       ylab = "Probability of detection",
       xlab = "Signal-to-noise ratio")
  
  ## Determine how many calls were simulated and extract distances to detectors
  n_call <- round(A * sum(D))
  cat(paste0(n_call, " bowhead whale calls were simulated.\n"))
  
  call_sample_index <- sample(1:n_grid, size = n_call, prob = D, 
                              replace = TRUE)
  distances <- distances[call_sample_index, ]
  calls <- cov_density[call_sample_index, ]
  
  ## Simulate sound characteristics for every call
  source_levels <- rtruncnorm(n_call,
                              a = par_sl["lower"],
                              b = par_sl["upper"],
                              mean = par_sl["mu"],
                              sd = par_sl["sd"])
  
  # Simulate mean noise level for every call
  noise <- rtruncnorm(n_call,
                      a = par_noise["lower"],
                      mean = par_noise["mu"],
                      sd = par_noise["sd"])
  # For every noise level, create slight deviations for all detectors
  noise <- matrix(rep(noise, each = n_det), nrow = n_call, ncol = n_det,
                  byrow = TRUE)
  noise <- noise + matrix(rnorm(length(noise), 0, 2), nrow = n_call,
                          ncol = n_det, byrow = TRUE)

  # For every call, derive received level and add measurement error
  received_levels <- matrix(source_levels, nrow = n_call, ncol = n_det,
                            byrow = FALSE) -
    par_rl["beta"] * log(distances, base = 10) + # transmission loss
    matrix(rnorm(n_call * n_det, 0, par_rl["sigma"]), nrow = n_call) # measurement error

  # Derive the signal to noise ratio
  snr <- received_levels - noise
  snr[snr < 0] <- 0
  
  ######################## Create detection histories ##########################
  
  ## Derive detection probabilities based on snr (set all negative snr to zero

  if (det_function == "half-normal") {
    det_probs <- .gHN(distances = distances, par = par_det)
  } else {
    det_probs <- .gSNR(snr = snr, par = par_det, type = det_function)
  }
  
  # Create detection history consisting of 1's and 0's
  det_hist <- det_probs > runif(n_call * n_det, min = 0, max = 1)
  mode(det_hist) <- "numeric"
  
  # Create index for whether min_no_detections was satisfied
  ENOUGH_DETS <- rowSums(det_hist) >= min_no_detections
  detected_calls <- calls[ENOUGH_DETS, ]
  
  # ## PLOTTING FOR TESTING
  # plot(cov_density[, c("x", "y")])
  # points(detectors, col = "red", pch = 16)
  
  # Print result of simulation
  cat(paste0(sum(ENOUGH_DETS), " bowhead whale calls were detected at least ", 
             min_no_detections, " times.\n"))    
  
  ####################### Create bearing histories #############################
  # Derive bearings for all calls with enough detections
  bearings <- apply(detectors[, c("x", "y")], 1, function(det) {
    x <- detected_calls$x - det[1]
    y <- detected_calls$y - det[2]
    return(coord2rad(x, y, control.circular = list(template = "geographics", 
                                                   units = "degrees",
                                                   modulo = "2pi")))
  })

  # ## PLOTTING FOR TESTING
  # plot(cov_density[, c("x", "y")])
  # points(detectors, col = "red", pch = 16)
  # points(detected_calls[2, c("x", "y")], col = "blue", pch = 16)
  
  # Add the measurement error
  errors <- matrix(circular::rvonmises(n = length(bearings), 
                             mu = circular(0, template = "geographics", 
                                           modulo = "2pi", units = "degrees"), 
                             kappa = kappa),
                   ncol = n_det, nrow = nrow(bearings))
  bearings <- bearings + errors
  bearings <- circular(bearings, template = "geographics", 
                       modulo = "2pi", units = "degrees")

  ################# Add the data for detected calls to dat #####################
  
  ## Create output data
  det_hist <- det_hist[ENOUGH_DETS, ]
  NO_DETECTION <- !det_hist
  
  received_level_hist <- received_levels[ENOUGH_DETS, ]
  received_level_hist[NO_DETECTION] <- NA

  noise_call <- noise[ENOUGH_DETS, ]
  noise_random <- noise[sample(1:n_call, 1000, replace = TRUE), ]

  bearing_hist <- bearings
  bearing_hist[NO_DETECTION] <- NA
  
  # Only essential added right now, could output more info but not required.
  dat <- list(det_hist = det_hist,
              received_levels = received_level_hist,
              noise_call = noise_call,
              noise_random = noise_random,
              bearings = bearing_hist)
  
  return(dat)
}
  