# Run for test runs ============================================================
# par <- par_start
# rm(list = setdiff(ls(), c("dat", "par")))
# method = "L-BFGS-B"
# maxit = 100
# TRACE = TRUE
# LSE = TRUE
# ==============================================================================

#' Fits an acoustic spatial capture recapture model to bowhead whale call data 
#' collected by Greeneridge from 2007 to 2014. The method requires information
#' on the observed detection histories, received levels, and noise recordings,
#' and information on the latent spatial origins and source levels of the calls. 
#' 
#' Fits an acoustic spatial capture-recapture model to the data, using a list of 
#' start values. Functionality was added to allow for different optimising 
#' methods, different number of iterations, to trace to progress or not, and 
#' to use the likelihood with or without the LSE-approximation. It is not 
#' possible to provide start values for the parameters of regression splines, 
#' but all other parameters require start values.
#' 
#' @param dat a list required and correctly named data
#' @param par a list of parameters to be estimated 
#' @param method which optimisation method to use (see \code{optim()} for more 
#' information)
#' @param maxit the maximum number of iterations (defaults to \code{100})
#' @param TRACE whether it should print the parameters values that are tried 
#' (defaults to \code{TRUE})
#' @param LSE whether it should use the log-sum-exp approximations (defaults to
#' \code{TRUE})
#'
#' @return An \code{S3.Object} of class 'bwASCR_model' containing the results of 
#' the fitting, the raw results from \code{optim()} and the original data inputs 
#' for completeness.
#'
#' @examples
#' bwASCR(...)
#'
#' @export
bwASCR <- function(dat, par, method = "L-BFGS-B", maxit = 100, TRACE = TRUE,
                   LSE = TRUE) {
  
  # REMEMBER TO UPDATE CERTAIN BITS DEPENDING ON WHICH llk...R SCRIPT WILL BE USED
  # LINES TO UPDATE: 80 (which pars to include), 159 (which confidence bounds to use),
  #                  146 (which llk....R script to use).
  
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
  
  COMPLETE_DATA <- all(names(data) %in% c("det_hist",
                                          "detectors",
                                          "cov_density",
                                          "received_levels",
                                          # "noise_call",
                                          # "noise_random",
                                          # "source_levels",
                                          "A_x",
                                          # "A_s",
                                          "f_density",
                                          "det_function",
                                          # "bearings",
                                          "min_no_detections"))
  
  if (!COMPLETE_DATA) {stop("'dat' argument is incomplete!")}
  
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
  
  # Set mixture for bearings to false if bearings are not used to avoid NULL comparison in if statements
  if (!USE_BEARINGS) {dat$BEAR_MIXTURE <- FALSE}
  BEAR_MIXTURE <- dat$BEAR_MIXTURE
  
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
  
  if (dat$SINGLE_SL) {
    par <- par[!names(par) %in% c("sd_s")] # USE ONLY FOR SINGLE SOURCE LEVEL
  }
  
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
  
  # Not sure if scaling is necessary ===========================================
  cov_density_scaled <- subset(dat$cov_density, select = -area)
  # cov_density_scaled[, -c(1, 2)] <- scale(subset(cov_density_scaled, 
  #                                                select = -c(long, lat)))
  # ============================================================================
  
  # Turn the bearings to radians, as the log likelihood function uses radians 
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
    fn <- .llkParallelSmooth
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
  
  # ============================================================================
  # Define bounds for optim()
  # ============================================================================
  if (dat$SINGLE_SL & dat$det_function == "simple" & USE_BEARINGS & !BEAR_MIXTURE) {
    lower_bounds <- c(g0 = -10, 
                      kappa = log(0), # -Inf
                      beta_r = log(1), 
                      sd_r = log(0.01),
                      mu_s = log(50), 
                      rep(-Inf, 100)) # this gives 0 on log and logit link
    upper_bounds <- c(g0 = 10, 
                      kappa = log(10000),
                      beta_r = log(50), 
                      sd_r = log(50),
                      mu_s  = log(300), 
                      rep(Inf, 100)) # this gives 2.7e10 on log and 1 on logit
  } else if (!dat$SINGLE_SL & dat$det_function == "simple" & USE_BEARINGS & !BEAR_MIXTURE) {
    lower_bounds <- c(g0 = -10, 
                      kappa = log(0), # -Inf
                      beta_r = log(1), 
                      sd_r = log(0.01),
                      mu_s = log(50), 
                      sd_s = log(0.01),
                      rep(-Inf, 100)) # this gives 0 on log and logit link
    upper_bounds <- c(g0 = 10, 
                      kappa = log(10000),
                      beta_r = log(50), 
                      sd_r = log(50),
                      mu_s  = log(300), 
                      sd_s = log(50),
                      rep(Inf, 100)) # this gives 2.7e10 on log and 1 on logit
  } else if (dat$SINGLE_SL & dat$det_function == "simple" & USE_BEARINGS & BEAR_MIXTURE) {
    lower_bounds <- c(g0 = -10, 
                      kappa_low = log(0), # -Inf
                      kappa_high = log(0), # -Inf
                      mix_par = -10,
                      beta_r = log(1), 
                      sd_r = log(0.01),
                      mu_s = log(50), 
                      rep(-Inf, 100)) # this gives 0 on log and logit link
    upper_bounds <- c(g0 = 10, 
                      kappa_low = log(10000),
                      kappa_high = log(10000),
                      mix_par = 10,
                      beta_r = log(50), 
                      sd_r = log(50),
                      mu_s  = log(300), 
                      rep(Inf, 100)) # this gives 2.7e10 on log and 1 on logit
  } else if (!dat$SINGLE_SL & dat$det_function == "simple" & USE_BEARINGS & BEAR_MIXTURE) {
    lower_bounds <- c(g0 = -10, 
                      kappa_low = log(0), # -Inf
                      kappa_high = log(0), # -Inf
                      mix_par = -10,
                      beta_r = log(1), 
                      sd_r = log(0.01),
                      mu_s = log(50), 
                      sd_s = log(0.01),
                      rep(-Inf, 100)) # this gives 0 on log and logit link
    upper_bounds <- c(g0 = 10, 
                      kappa_low = log(10000),
                      kappa_high = log(10000),
                      mix_par = 10,
                      beta_r = log(50), 
                      sd_r = log(50),
                      mu_s  = log(300), 
                      sd_s = log(50),
                      rep(Inf, 100)) # this gives 2.7e10 on log and 1 on logit
  } else if (dat$SINGLE_SL & dat$det_function == "simple" & !USE_BEARINGS) {
    lower_bounds <- c(g0 = -10, 
                      beta_r = log(1), 
                      sd_r = log(0.01),
                      mu_s = log(50), 
                      rep(-Inf, 100)) # this gives 0 on log and logit link
    upper_bounds <- c(g0 = 10, 
                      beta_r = log(50), 
                      sd_r = log(50),
                      mu_s  = log(300), 
                      rep(Inf, 100)) # this gives 2.7e10 on log and 1 on logit
  } else if (!dat$SINGLE_SL & dat$det_function == "simple" & !USE_BEARINGS) {
    lower_bounds <- c(g0 = -10, 
                      beta_r = log(1), 
                      sd_r = log(0.01),
                      mu_s = log(50), 
                      sd_s = log(0.01),
                      rep(-Inf, 100)) # this gives 0 on log and logit link
    upper_bounds <- c(g0 = 10, 
                      beta_r = log(50), 
                      sd_r = log(50),
                      mu_s  = log(300), 
                      sd_s = log(50),
                      rep(Inf, 100)) # this gives 2.7e10 on log and 1 on logit
  }
  
  
  ######################## Fit using optim() ###################################
  
  fit_duration <- system.time({
    # Use bounds to limit search space to sensible space
    result <- optim(par = par, 
                    fn = fn, 
                    method = method, 
                    hessian = USE_MLE_SD,
                    control = list(maxit = maxit, 
                                   fnscale = -1, 
                                   trace = TRACE,
                                   factr = 1e10, # The convergence factor (lower is more precise). defaults to 1e7, which is a tolerance of roughly 1e-8
                                   REPORT = 1),
                    lower = lower_bounds, 
                    upper = upper_bounds,
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
    N <- .totalN(D = D$density, A = dat$A_x$area, n = n,
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
  
  if (TRACE) {
    output <- list(N = N, estimates = estimates, real = real, 
                   n_par = K, AIC = aic, AICc = aicc, BIC = bic,
                   start_values = start_values, method = method, hessian = hess,
                   covariance_matrix = covariance_matrix, density = D, 
                   det_function = dat$det_function, design_matrix = design, 
                   ESA = esa, optim_result = result, fit_duration = fit_duration, 
                   data = dat, par_hist = read.csv(paste0("parameter_history_", 
                                                          dat$f_density[3], ".csv"), 
                                                   header = TRUE))
  } else {
    output <- list(N = N, estimates = estimates, real = real, 
                   n_par = K, AIC = aic, AICc = aicc, BIC = bic,
                   start_values = start_values, method = method, hessian = hess,
                   covariance_matrix = covariance_matrix, density = D, 
                   det_function = dat$det_function, design_matrix = design, 
                   ESA = esa, optim_result = result, fit_duration = fit_duration, 
                   data = dat, par_hist = NA)
  }
  
  class(output) <- "bwASCR_model"
  
  return(output)
}