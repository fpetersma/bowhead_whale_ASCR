# # for test runs
# par <- par_start
# rm(list = setdiff(ls(), c("dat", "par")))
# method = "L-BFGS-B"
# maxit = 100
# TRACE = TRUE
# LSE = TRUE

bwASCR <- function(dat, par, method = "L-BFGS-B", maxit = 100, TRACE = TRUE,
                   LSE = TRUE) {
  
  # Description: 
  #   An acoustic spatial capture-recapture model to the data, using a list of 
  #   start values. Functionality was added to allow for different optimising 
  #   methods, different number of iterations, to trace to progress or not, and 
  #   to use the likelihood with or without the LSE-trick. It is not possible to 
  #   provide start values for the parameters of regression splines, but all 
  #   other parameters require start values.
  
  # REMEMBER TO UPDATE CERTAIN BITS DEPENDING ON WHICH llk...R SCRIPT WILL BE USED
  # LINES TO UPDATE: 80 (which pars to include), 159 (which confidence bounds to use),
  #                  146 (which llk....R script to use).
  
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
  library(raster)
  
  # Input checks
  if (class(dat) != "list") {
    stop()
  }
  
  # Use normal MLE standard error?
  USE_MLE_SD <- FALSE
  
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
  if (dat$det_function == "janoschek" | dat$det_function == "logit" | dat$det_function == "probit") {
    CORRECT_PAR <- all(c("U", "B", "Q") %in% names(par[["par_det"]]))
    if (!CORRECT_PAR) {
      stop("Incorrect start parameters specified for the SNR detection function.")
    }
  } else if (dat$det_function == "half-normal") {
    CORRECT_PAR <- all(c("g0", "sigma") %in% names(par[["par_det"]]))
  } else if (dat$det_function == "simple") {
    CORRECT_PAR <- all(c("g0") %in% names(par[["par_det"]]))
  } else {
    stop(paste0("Detection function specification is ", dat$det_function,
                ", but should be 'simple', 'janoscheck', 'logit', 'probit'  or 'half-normal'"))
  }
  
  # Turn par into a named vector with the correct names (optim() requires a vector)
  par <- unlist(par)
  names(par) <- gsub(".*\\.", "", names(par))
  ## TO KEEP SOME PARAMETERS FIXED, USE LINE BELOW TO SELECT WHICH PARAMETERS TO ESTIMATE
  # par <- par[names(par) %in% c("g0", "sigma" , "kappa", "(Intercept)", "dist_to_coast", "dist_to_coast2")] # USE ONLY FOR SIMULATIONS
  # par <- par[!names(par) %in% c("sd_r")] # USE ONLY FOR SINGLE SOURCE LEVEL
  
  # Create a matrix of distances
  distances <- pointDistance(p1 = dat$cov_density[, c("long", "lat")], 
                             p2 = dat$detectors[, c("long", "lat")], 
                             lonlat = TRUE)
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
  
  # Not sure if scaling is necessary
  cov_density_scaled <- subset(dat$cov_density, select = -area)
  cov_density_scaled[, -c(1, 2)] <- scale(subset(cov_density_scaled, 
                                                 select = -c(long, lat)))
  
  # Turn the bearings to radians, as the llk calculation uses radians 
  if (USE_BEARINGS) {
    bearings_deg <- circular(dat$bearings,
                             units = "degrees",
                             template = "geographics") # 'geographics' means clockwise rotation with 0 at north
    bearings_rad <- conversion.circular(bearings_deg)
    dat[["bearings_rad"]] <- bearings_rad # Add bearings as circular in radians to dat
    
    # Get degrees bearings from detectors to all grid points 
    grid_bearings <- t(apply(dat$cov_density[, c("long", "lat")], 1, function(grid_point) {
      bear <- geosphere::bearing(p1 = as.matrix(dat$detectors[, c("long", "lat")]),
                                 p2 = grid_point)
    }))
    # First convert to 'circular' 
    grid_bearings <- circular(grid_bearings, template = "geographics", 
                              modulo = "2pi", units = "degrees")
    # Then convert to radians
    grid_bearings <- conversion.circular(grid_bearings, units = "radians")
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
    
    # par <- c(par, gam_par)
    par <- c(par, gam_fit$coefficients[names(gam_fit$coefficients) != "(Intercept)"])
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
  
  if (TRACE) {
    write.table(matrix(par, nrow = 1, byrow = TRUE), 
                file = paste0("parameter_history_", dat$f_density[3], ".csv"),
                sep = ",",
                row.names = FALSE,
                col.names = names(par))
  }
  
  ######################## Fit using optim() ###################################
  
  fit_duration <- system.time({
    # Use bounds to limit search space to sensible space
    result <- optim(par = par, fn = fn, method = method, hessian = USE_MLE_SD,
                    control = list(maxit = maxit, fnscale = -1, trace = TRACE, 
                                   REPORT = 1, factr = 1e10),
                    # lower = c(U = -5, B = -5, Q = log(0.01), kappa = log(5), 
                    #           beta_r = log(1), mu_s = log(50), sd_r = log(0.01),
                    #           sd_s = log(0.01), rep(-Inf, 100)), # this gives 0 on log and logit link
                    # upper = c(U = 5, B = 5, Q = log(10000), kappa = log(100), 
                    #           beta_r = log(50), mu_s  = log(300), sd_r = log(50),
                    #           sd_s = log(50), rep(Inf, 100)), # this gives 2.7e10 on log and 1 on logit
                    dat = dat)
  })
  
  ############# Create additional output with estimated parameters #############
  # Calculating the final density
  par_density <- result$par[colnames(dat$design)]
  D <- cbind(dat$cov_density[, c("long", "lat")], 
             density = .densityGAM(dat$design, dat$gam_fit, par_density))
  
  if (USE_MLE_SD) {
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
  } else {
    hess <- NA
    covariance_matrix <- NA
    
    K <- length(par) # number of parameters
    n <- nrow(dat$det_hist) # number of calls detected
    
    # Derive the estimated pars
    link <- rep("log", K)
    link[1] <- "logit"
    if (dat$det_function == "janoschek") {
      link[names(par) == "Q"] <- "log*"
    }
    pars <- result$par
    estimates <- data.frame(link = link, estimate = round(pars, 6))
    
    # Derive the pars on the real scale
    real_pars <- vector(mode = "numeric")
    for (i in seq_along(pars)) {
      if (link[i] == "logit") {
        real_pars[i] <- 1 / (1 + exp(-pars[i]))
      } else if (link[i] == "log*") {
        real_pars[i] <- exp(pars[i]) + 1
      } else {
        real_pars[i] <- exp(pars[i])
      }
    }
    real <- data.frame(link = link, estimate = round(real_pars, 6))
    row.names(real) <- row.names(estimates)
  }
  
  design <- model.matrix(gam_fit)
  
  if (USE_MLE_SD) {
    N <- .total_N(D = D, A = dat$A_x$area, n = n,
                  covariance_matrix = covariance_matrix, design = design)
  } else {
    N <- c(estimate = sum(dat$A_x$area * D$density))
  }
  
  aic <- -2 * result$value + 2 * K # Burnham & Anderson (2002)
  aicc <- aic + 2 * K * (K + 1) / (n - K - 1) # Burnham & Anderson (2002)
  bic <- -2 * result$value + K * log(n)
  
  # ESA not calculated here anymore
  esa <- NA
  
  ## Finalise the process
  cat("Finished the SCR fitting.\n")
  cat(paste0("Total number of calls was estimated at ", round(N["estimate"]), "\n"))
  
  output <- list(N = N, estimates = estimates, real = real, 
                 n_par = K, AIC = aic, AICc = aicc, BIC = bic,
                 start_values = start_values, method = method, hessian = hess,
                 covariance_matrix = covariance_matrix, density = D, 
                 det_function = dat$det_function, design_matrix = design, 
                 ESA = esa, optim_result = result, fit_duration = fit_duration, 
                 data = dat, par_hist = read.csv("parameter_history.csv", 
                                                 header = TRUE))
  class(output) <- "bwASCR_model"
  
  return(output)
}