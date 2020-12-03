# # for test runs
# par <- par_start
# rm(list = setdiff(ls(), c("dat", "par")))
# method = "L-BFGS"
# maxit = 100
# TRACE = TRUE
# LSE = TRUE

bwASCR <- function(dat, par, method = "L-BFGS", maxit = 100, TRACE = TRUE,
                   LSE = TRUE) {
  
  # Description: 
  #   n acoustic spatial capture-recapture model to the data, using a list of 
  #   start values. Functionality was added to allow for different optimising 
  #   methods, different number of iterations, to trace to progress or not, and 
  #   to use the likelihood with or without the LSE-trick. It is not possible to 
  #   provide start values for the parameters of regression splines, but all 
  #   other parameters require start values.
  #   
  
  # Inputs:
  #   dat       - [list] 
  #   param     - [list] 
  #   method    - [character] 
  #   maxit     - [scalar] 
  #   trace     - [logical] 
  #   LSE       - [logical]
  
  # Outputs:
  #   output    - [S3.Object] containing the results of the fitting, the raw 
  #               results from optim() and the original data inputs for 
  #               completeness.
  
  ################## Load libraries and perform input checks ###################
  library(dplyr)
  library(mgcv)
  library(matrixStats)
  library(circular)
  
  source("Scripts/Bowhead Whales/HTLikeEstimator.R")
  
  # Input checks
  if (class(dat) != "list") {
    stop()
  }
  ##############################################################################
  ########################## Start the fitting #################################
  ##############################################################################
  
  cat(paste0("Fitting the SCR model using the ", method, " method.\n"))
  
  ############ Perform all necessary checks and preparations ###################
  
  start_values <- par  # Save start values
  
  USE_BEARINGS <- "par_bear" %in% names(par) # Check whether bearings are provided
  USE_RL <- all(c("par_sl", "par_rl") %in% names(par)) # Check whether received levels are provided
  
  dat[["USE_BEARINGS"]] <- USE_BEARINGS
  dat[["USE_RL"]] <- USE_RL
  
  # Check which detection function is to be used, and make sure detection rate
  # or probability at distance = 0 is correctly named. 
  if (dat$det_function == "janoschek" | dat$det_function == "logistic") {
    CORRECT_PAR <- all(c("U", "B", "Q") %in% names(par[["par_det"]]))
    if (!CORRECT_PAR) {
      stop("Incorrect start parameters specified for the SNR detection function.")
    }
  } else if (dat$det_function == "half-normal") {
    CORRECT_PAR <- all(c("g0", "sigma") %in% names(par[["par_det"]]))
  } else {
    stop(paste0("Detection function specification is ", dat$det_function,
                ", but should be 'janoscheck', 'logistic' or 'half-normal'"))
  }
  
  # Turn par into a named vector with the correct names (optim() requires a vector)
  par <- unlist(par)
  names(par) <- gsub(".*\\.", "", names(par))
  ## TO KEEP SOME PARAMETERS FIXED, USE LINE BELOW TO SELECT WHICH PARAMETERS TO ESTIMATE
  # par <- par[names(par) %in% c("g0", "sigma" , "kappa", "(Intercept)", "dist_to_coast", "dist_to_coast2")] # USE ONLY FOR SIMULATIONS
  # par <- par[names(par) %in% c("U", "B", "Q", "mu_s", "beta_r", "(Intercept)")] # USE ONLY FOR SIMULATIONS
  
  # Create a matrix of distances
  distances <- apply(dat$detectors, 1, function(z) {
    y <- t(dat$cov_density[, c("x", "y")]) - z
    y <- y ^ 2
    return(sqrt(colSums(y)))
  })
  dat[["distances"]] <- distances  # Add matrix of distances to dat
  
  # Extract variables from f_density, remove 'D' and add '(Intercept)'
  vars <- all.vars(dat$f_density)
  vars <- c("(Intercept)", vars[-1])
  
  # Check whether density is constant
  if (length(vars) == 1) {
    CONSTANT_DENSITY <- FALSE
  } else {
    CONSTANT_DENSITY <- FALSE
  }
  dat["CONSTANT_DENSITY"] <- CONSTANT_DENSITY
  
  # not sure if scaling is necessary
  cov_density_scaled <- dat$cov_density
  cov_density_scaled[, -c(1, 2)] <- scale(select(dat$cov_density, -c("x", "y")))
  
  if (USE_BEARINGS) {
    bearings_deg <- circular(dat$bearings, 
                             units = "degrees", 
                             template = "geographics") # 'geographics' means clockwise rotation with 0 at north
    bearings_rad <- conversion.circular(bearings_deg)
    dat[["bearings_rad"]] <- bearings_rad # Add bearings as circular in radians to dat
    
    # Add bearings of all grid points to all detectors to dat
    grid_bearings <- apply(dat$detectors[, c("x", "y")], 1, function(det) {
      x <- dat$cov_density$x - det[1]
      y <- dat$cov_density$y - det[2]
      return(coord2rad(x, y, control.circular = list(template = "geographics")))
    })
    dat[["grid_bearings"]] <- grid_bearings
  }
  
  # Create a gam object (not necessary if CONSTANT_DENSITY == FALSE?)
  gam_fit <- gam(dat$f_density, 
                 data = cbind(D = 0, cov_density_scaled))
  if (length(gam_fit$smooth) != 0) {
    k <- gam_fit$smooth[[1]]$bs.dim
    gam_par <- rep(0, k - 1)
    dat[["smooth_terms"]] <- paste0(gam_fit$smooth[[1]]$label, ".", 1:(k - 1))
    names(gam_par) <- dat$smooth_terms
    
    par <- c(par, gam_par)
  } else {dat[["smooth_terms"]] <- NULL}
  
  dat[["design"]] <- model.matrix(gam_fit)
  dat[["gam_fit"]] <- gam_fit
  
  # Add whether trace is TRUE/FALSE
  dat[["TRACE"]] <- TRACE
  
  if (LSE) {
    fn <- .llkLSEParallelSmooth
  } else {
    stop("Non LSE version not available yet.")
    #fn <- .loglikelihood_with_bearings
  }
  
  # Remove density from par if density is constant STILL NEEDS IMPLEMENTATION
  if (CONSTANT_DENSITY) {
    par <- par[!(names(par) %in% vars)]
  }
  
  ######################## Fit using optim() ###################################
  
  # Use bounds to limit search space to sensible space
  result <- optim(par = par, fn = fn, method = method, hessian = TRUE,
                  control = list(maxit = maxit, fnscale = -1, trace = TRACE, 
                                 REPORT = 1, factr = 1e11),
                  # lower = c(U = -5, B = -5, Q = log(1), kappa = log(1), beta_r = log(10), sd_r = log(0.1), mu_s = log(70),  sd_s = log(1)), # this gives 0 on log and logit link
                  # upper = c(U = 5, B = 5, Q = log(10), kappa = log(100), beta_r = log(20), sd_r = log(10),   mu_s  =log(130) , sd_s = log(20) ), # this gives 2.7e10 on log and 1 on logit
                  # # lower = c(g0 = -5, sigma = log(5000), kappa = log(1)),#, beta_r = log(10), sd_r = log(0.1), mu_s = log(70),  sd_s = log(1)), # this gives 0 on log and logit link
                  # upper = c(g0 = 5, sigma = log(50000), kappa = log(100)),#, beta_r = log(20), sd_r = log(10),   mu_s  =log(130) , sd_s = log(20) ), # this gives 2.7e10 on log and 1 on logit
                  dat = dat)
  # return(result)
  
  if (CONSTANT_DENSITY) {
    ## Deriving the Horvitz-Thompson-like estimator for constant density
    n_call <- nrow(dat$det_hist)
    if (dat$det_function == "janoschek" | dat$det_function == "logistic") {
      par_det <- c(exp(result$par["U"]) / (1 + exp(result$par["U"])),
                   exp(result$par["B"]),
                   exp(result$par["Q"]) + 1)
      par_sl <- exp(result$par[c("mu_s", "sd_s")])
      beta_r <- exp(result$par["beta_r"])
      
      d <- HTLikeEstimator(distances = distances, n_call = n_call, A = dat$A,
                           par_det = par_det, det_function = dat$det_function, 
                           min_no_detections = dat$min_no_detections, 
                           noise = dat$noise_call, beta_r = beta_r, 
                           source_levels = dat$source_levels, par_sl = par_sl)
    } else if (dat$det_function == "half-normal") {
      par_det <- c(exp(result$par["g0"]) / (1 + exp(result$par["g0"])),
                   exp(result$par["sigma"]))
      d <- HTLikeEstimator(distances = distances, n_call = n_call, A = dat$A,
                           par_det = par_det, det_function = dat$det_function, 
                           min_no_detections = dat$min_no_detections)
    }
    return(list(result = result, constant_density = d))
    
  } else {
    return(list(result = result))
  }
  
  ############# Create additional output with estimated parameters #############
  # Calculating the final density
  par_density <- result$par[colnames(dat$design)]
  D <- cbind(dat$cov_density[, c("x", "y")], 
             density = .densityGAM(dat$design, dat$gam_fit, par_density))
  
  hess <- result$hessian
  
  # When maximising a likelihood then the covariance matrix of the 
  # estimates is (asymptotically) the inverse of information matrix.
  K <- length(par) # number of parameters
  n <- nrow(dat$det_hist) # number of calls detected
  covariance_matrix <- NULL
  estimates <- NULL
  
  # Question is, since we are optimising, do we want the negative of the Hessian
  # or is the hessian itself the information matrix (since fn is already 
  # multiplied by -1)?
  try({
    covariance_matrix <- solve(-hess)
    se <- sqrt(diag(covariance_matrix))  
    
    pars <- result$par
    
    lower <- pars - 1.96 * se
    upper <- pars + 1.96 * se
    
    link <- rep("log", K)
    link[1] <- "logit"
    if (dat$det_function == "janoschek") {
      link[names(par) == "Q"] <- "log*"
    }
    # if (dat$det_function %in% c("HN", "HR")) {
    #   link[names(par) == "g0"] <- "logit"
    # }
    # if (USE_SNR) {
    #   link[names(param) == "b1"] <- "identity"
    # }
    # link <- c("log", "logit", rep("log", length(result$par)-2))
    
    estimates <- data.frame(link = link, estimate = round(pars, 6), 
                            se = round(se, 6), lower = round(lower, 6), 
                            upper = round(upper, 6))
    
    # Create df for real parameter estimates evaluated at base values (i.e. 0)
    real_pars <- vector(mode = "numeric")
    real_lower <- vector(mode = "numeric")
    real_upper <- vector(mode = "numeric")
    for (i in seq_along(pars)) {
      if (link[i] == "logit") {
        real_pars[i] <- 1 / (1 + exp(-pars[i]))
        real_lower[i] <- 1 / (1 + exp(-lower[i]))
        real_upper[i] <- 1 / (1 + exp(-upper[i]))
      } else if (link[i] == "log*") {
        real_pars[i] <- exp(pars[i]) + 1
        real_lower[i] <- exp(lower[i]) + 1
        real_upper[i] <- exp(upper[i]) + 1
      } else {
        real_pars[i] <- exp(pars[i])
        real_lower[i] <- exp(lower[i])
        real_upper[i] <- exp(upper[i])
      }
    }
    real <- data.frame(link = link, estimate = round(real_pars, 6), 
                       lower = round(real_lower, 6), 
                       upper = round(real_upper, 6))
    row.names(real) <- row.names(estimates)
  })
  
  design <- model.matrix(gam_fit)
  
  N <- .total_N(D = D, A = dat$A, n = n,
                covariance_matrix = covariance_matrix, design = design)
  
  aic <- -2 * result$value + 2 * K # Burnham & Anderson (2002)
  aicc <- aic + 2 * K * (K + 1) / (n - K - 1) # Burnham & Anderson (2002)
  bic <- -2 * result$value + K * log(n)
  
  # # Calculate effective sample area (a)
  # if (det_function == "HHN") {
  #   det_probs <- .p(distances = distances, det0 = estimates["lambda0", 2], 
  #                   sigma = estimates["sigma", 2], det_function = det_function)
  # } else if (det_function == "HN") {
  #   det_probs <- .p(distances = distances, det0 = estimates["g0", 2], 
  #                   sigma = estimates["sigma", 2], det_function = det_function)
  # }
  # p. <- .detected(det_probs = det_probs, min_no_detections = min_no_detections)
  # esa <- sum(p. * D$density) / sum(D$density)
  #
  # rm(det_probs, p.)
  
  
  esa <- "STILL NEEDS WORK"
  
  ## Finalise the process
  cat("Finished the SCR fitting.\n")
  cat(paste0("Total number of calls was estimated at ", round(N["estimate"]), "\n"))
  
  
  output <- list(N = N, estimates = estimates, real = real, 
                 n_par = K, AIC = aic, AICc = aicc, BIC = bic,
                 start_values = start_values, method = method, hessian = hess,
                 covariance_matrix = covariance_matrix, density = D, 
                 det_function = dat$det_function, design_matrix = design, 
                 ESA = esa, optim_result = result, data = dat)
  class(output) <- "bwASCR_model"
  
  return(output)
}