HTLikeEstimator <- function(distances, n_call, par_det, det_function, A = A,
                            min_no_detections, noise = NULL, beta_r = NULL,
                            source_levels = NULL, par_sl = NULL) {
  # Derive the Horvitz-Thompson-like estimator for density, when density is 
  # assumed to be homogeneous in space.
  
  # Input: 
  #   distances is the covariate for the detection fucntion
  #   n is the number of detection animals
  #   par_det is a named vector containing the detection function parameters
  #   det_function is the name of the detection function
  #   noise is a matrix with a noise vector for every detected call (only 
  #   required when det_function == "jano")
  if (det_function == "half-normal") {
    # Derive the detection probabilities based on distances
    det_probs <- .gHN(distances = distances, par = par_det)
    
    # Derive probabilities that a call was detected at least 
    # min_no_detections times for all source levels
    p. <- .detected(probs = det_probs, min_no_detections = min_no_detections)
    
    # Sum all the p. values and multiply by A to get the ESA 
    ESA <- sum(p.) * A
    
    # Horvitz-Thompson-like estimator D = n / ESA
    D <- n_call / ESA
    
    return(D)
    
  } else if (det_function == "jano") {
    if (is.null(noise)) {
      stop("Noise data is required to derive the Horvitz-Thompson-like estimator.")
    }
    n_sl <- length(source_levels)
    n_det <- ncol(noise)
    n_grid <- nrow(distances)
    
    mu_s <- par_sl["mu_s"]
    sd_s <- par_sl["sd_s"]
    
    ESA <- apply(t(noise), 2, function(c) {
      # Get the expected SNR for all noise samples, grid points and detectors
      E_snr <- array(data = NA,
                     dim = c(n_sl, n_grid, n_det),
                     dimnames = list("source_level" = source_levels, 
                                     "grid_point" = 1:n_grid, 
                                     "detector" = 1:n_det))
      
      # Create matrix with identical noise in every row
      c_m <- matrix(c, nrow = n_grid, ncol = n_det, byrow = TRUE)
    
      for (sl in source_levels) {
        sl_m <- matrix(sl, nrow = n_grid, ncol = n_det)
        # Derive the expected snr
        E_snr[as.character(sl), , ] <- sl_m - beta_r * 
          log(distances, base = 10) - c_m
      }
      
      # Set E_snr to 0 if negative
      E_snr[E_snr < 0] <- 0
    
      # Clean environment
      rm(sl_m, sl) 
      
      # Derive the associated detection probabilities 
      det_probs <- .gJano(snr = E_snr, par = par_det)

      # Derive probabilities that a call was detected at least 
      # min_no_detections times for all source levels
      p. <- t(apply(det_probs, 1, function(probs) {
        .detected(probs = probs, min_no_detections = min_no_detections)
      }))
      
      # Run some test whether the data is correct
      if (any(det_probs < 0)) {
        print("At least one element in det_probs is negative!")
      }
      if (any(p. < 0)) {
        print("At least one element in p. is negative!")
      }
      
      # Create a matrix of probability densities for source level
      probs_sl <- dnorm(source_levels, mean = mu_s, sd = sd_s)
      probs_sl_m <- matrix(probs_sl, nrow = n_sl, ncol = n_grid, 
                           byrow = FALSE)
      
      # Create a matrix of all products of D, p. and f(s)
      probs_per_grid_sl <- p. * probs_sl_m
      
      # Sum over all columns (grid points) 
      probs_per_sl <- Rfast::colsums(t(probs_per_grid_sl)) * A
      
      sum(probs_per_sl)
    })
    
    D <- sum(1 / ESA)
    
    return(D)
  }

}