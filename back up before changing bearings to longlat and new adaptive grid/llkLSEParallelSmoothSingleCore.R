#### Version of llkLSE_parallel.R, but now with a smooth fitted to the      ####
#### likelihoods for the source levels.

#### The big question is where the smooths should occur. My first intuition is 
#### to do it at the source level, since this is a nice univariate smooth.
#### At the noise level it would be a six dimensional smooth and for the grid
#### points is would be a two dimensional smooth. 

# # for test runs
# par <- par
# dat <- dat
# rm(list = setdiff(ls(), c("dat", "par")))

# The 'fast' loglikelihood (not really, but faster than otherwise)
.llkLSEParallelSmooth <- function(par, dat) {
  
  # Description:
  #   
  
  # Inputs:
  #   param       - [list] 
  #   dat         - [list] 
  
  # Outputs:
  #   total       - [scalar] the total log likelihood value.
  
  ###################### INPUT CHECKS STILL MISSING ############################
  require(Rfast)
  require(mgcv)
  library(parallel)
  library(snow)
  
  # Define very small error to be added to values
  cores <- 12
  error <- 0
  smooth <- "none" # other options are "loess" and "gam"
  
  ######################## Extract data from dat ###############################
  det_hist <- dat$det_hist
  cov_density <- dat$cov_density
  detectors <- dat$detectors
  A_x <- dat$A_x
  f_density <- dat$f_density
  det_function <- dat$det_function
  distances <- dat$distances
  gam_fit <- dat$gam_fit
  design <- dat$design
  min_no_detections <- dat$min_no_detections
  USE_BEARINGS <- dat$USE_BEARINGS
  USE_RL <- dat$USE_RL
  TRACE <- dat$TRACE
  CONSTANT_DENSITY <- dat$CONSTANT_DENSITY
  
  if (TRACE) {print(par)}
  
  if (USE_BEARINGS) {
    bearings_rad <- dat$bearings_rad
    grid_bearings <- dat$grid_bearings
  }
  if (USE_RL) {
    A_s <- dat$A_s
    source_levels <- dat$source_levels
    received_levels <- dat$received_levels
    noise_call <- dat$noise_call
    noise_random <- dat$noise_random
  }
  
  ################## Extract some important constants ##########################
  n_call <- nrow(det_hist) # number of calls
  n_det <- nrow(detectors) # number of detectors
  n_grid <- nrow(cov_density) # number of grid points
  if (USE_RL) {
    n_sl <- length(source_levels) # number of source levels
    n_noise <- nrow(noise_random) # number of random noise samples
  }
  
  if (!CONSTANT_DENSITY) {
    vars <- colnames(design) # parameters for density function
  }
  
  # Extract parameters and convert to the real scale (add errors to ensure 
  # correct domain)
  if (det_function == "janoschek" ) {
    U <- exp(par["U"]) / (1 + exp(par["U"])) # logit link 
    B <- exp(par["B"]) + error # log link
    Q <- exp(par["Q"]) + 1 #+ error # log link + 1
    par_det <- c(U, B, Q)
    names(par_det) <- c("U", "B", "Q") # make sure names are correct
  } else if (det_function == "logit" | det_function == "probit") {
    U <- exp(par["U"]) / (1 + exp(par["U"])) # logit link 
    B <- exp(par["B"]) + error # log link
    Q <- exp(par["Q"]) + error # log link 
    par_det <- c(U, B, Q)
    names(par_det) <- c("U", "B", "Q") # make sure names are correct
  } else if (det_function == "half-normal") {
    g0 <- exp(par["g0"]) / (1 + exp(par["g0"])) # logit link 
    sigma <- exp(par["sigma"]) # log link
    par_det <- c(g0, sigma)
    names(par_det) <- c("g0", "sigma") # make sure names are correct
  } else {stop("A PROBLEM!")}
  
  if (USE_BEARINGS) {
    kappa <- exp(par["kappa"]) # log link
  } 
  
  if (USE_RL) {
    beta_r <- exp(par["beta_r"]) # log link
    sd_r <- exp(par["sd_r"]) # log link
    
    mu_s <- exp(par["mu_s"]) # log link
    sd_s <- exp(par["sd_s"]) # log link
    lower_s <- 0
    upper_s <- Inf
  }
  
  if (smooth != "none" & USE_RL) {
    # Get a subset of source levels to fit a smooth to
    source_levels_all <- source_levels
    source_levels <- source_levels[seq(from = 1, to = 61, by = 6)]
    n_sl <- length(source_levels)
  }
  
  if (!CONSTANT_DENSITY) {
    # Check if all density parameters have starting values
    if (!all(vars %in% names(par))) {
      stop("Not all density parameters have starting values.")
    }
    par_dens <- par[vars]
  }
  ##############################################################################
  ###################### Start deriving the log likelihood #####################
  ##############################################################################
  
  ################ Derive probabilites for all combinations ####################
  
  if (!CONSTANT_DENSITY) {
    # Derive the density
    D <- .densityGAM(x = design, gam_fit = gam_fit, par = par_dens)
    if (any(D < 0)) {
      print("At least one element in D is negative!")
    }
  } else {
    D <- rep(1, n_grid)
  }
  
  if (USE_RL) {
    ################ Start parallelisation over source level #####################
    
    ### Creating a list of one dimension and than parallel lapply() ##############
    # Faster than foreach()
    # library(parallel)
    # library(snow)
    
    # system.time({
    if (!CONSTANT_DENSITY) {
      ##########################################################################
      ################ Derive the estimate for rate parameter ##################
      ##########################################################################
      ## Some experimenting showed that using 60 random noise samples, the
      ## parallel version took roughly 20secs, whereas a normal lapply took
      ## 60secs, roughly 3 times as long. This will only increase as the size
      ## or noise_random increases.
      
      ## Initiate parallel process leaving one core free
      # no_cores <- cores#detectCores() - 1
      # cl <- makeCluster(no_cores)
      
      # Export required data and functions to all clusters
      # hidden_functions <- c(".gSNR", ".detected", ".densityGAM", ".gHN")
      # clusterExport(cl, list = c(ls(), hidden_functions), envir = environment()) 
      
      # Turn noise random in a t(data frame) and then into a list to create a 
      # list with the noise for six detectors as a vector in every element.
      noise_random_list <- as.list(as.data.frame(t(noise_random)))
      
      # Start parallel process over noise_random
      # rates_noise <- parLapply(cl, noise_random_list, function(c) {
      rates_noise <- lapply(noise_random_list, function(c) {
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
        
        # # Set E_snr to 0 if negative
        # E_snr[E_snr < 0] <- 0
        # Clean environment
        rm(sl_m, sl) 
        
        # Derive the associated detection probabilities 
        if (det_function == "janoschek" | det_function == "logit" | det_function == "probit") {
          det_probs <- .gSNR(snr = E_snr, 
                             par = par_det,
                             type = det_function)
        } else if (det_function == "half-normal") {
          # Create same distance array for every source level
          distance_array <- array(data = NA,
                                  dim = c(n_sl, n_grid, n_det),
                                  dimnames = list("source_level" = source_levels, 
                                                  "grid_point" = 1:n_grid, 
                                                  "detector" = 1:n_det))
          for (i in 1:n_sl) {
            distance_array[i, , ] <- distances
          }
          det_probs <- .gHN(distances = distance_array, 
                            par = par_det)
        }
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
        
        # Create a matrix of density
        D_m <- matrix(D, nrow = n_sl, ncol = n_grid, byrow = TRUE)
        
        # Create a matrix of probability densities for source level
        probs_sl <- truncnorm::dtruncnorm(x = source_levels, a = lower_s, 
                                          b = upper_s, mean = mu_s, sd = sd_s)
        probs_sl_m <- matrix(probs_sl, nrow = n_sl, ncol = n_grid, 
                             byrow = FALSE)
        
        # Create a matrix of all products of D, p. and f(s)
        probs_per_grid_sl <- p. * D_m * probs_sl_m
        
        # Sum over all columns (grid points) 
        probs_per_sl <- Rfast::colsums(t(probs_per_grid_sl))
        
        if (smooth == "gam") {
          # Fit a smooth 
          smooth_fit <- mgcv::gam(probs_per_sl ~ s(source_levels)) # gam is more robust but slower
          pred <- mgcv::predict.gam(smooth_fit,
                                    newdata = data.frame(source_levels = source_levels_all))
        } else if (smooth == "loess") {
          smooth_fit <- loess(probs_per_sl ~ source_levels) # loess has edge effects
          pred <- predict(smooth_fit,
                          newdata = data.frame(source_levels = source_levels_all))
        } else if (smooth == "none") {pred <- probs_per_sl}
        
        return(sum(pred))
      }) 
      # Stop parallel process
      # stopCluster(cl)
      
      # Sum the output together and multiply by A to get the rate parameter 
      rate <- A_s * A_x * mean(unlist(rates_noise))
      
      # Derive the probability density for the number of detections given the rate
      llk_n <- dpois(n_call, rate, log = TRUE)
      if (llk_n == -Inf) llk_n <- -6000 # set to -6000 if -Inf since L-BFGS-B does
      # not accept non-real values
      
      # Clear environment
      rm("rate", "cl", "rates_noise", "no_cores", "hidden_functions", 
         "noise_random_list")
    } else {llk_n <- 0}
    # })
    ##############################################################################
    ############### Derive the conditional likelihood part #######################
    ##############################################################################
    
    ## Initiate parallel process leaving one core free
    # no_cores <- cores# detectCores() - 1
    # cl <- makeCluster(no_cores)
    
    # Export required data and functions to all clusters
    # hidden_functions <- c(".gSNR", ".detected", ".densityGAM", ".gHN")
    # clusterExport(cl, list = c(ls(), hidden_functions), envir = environment()) 
    
    # Start parallel process over number of calls 
    # cond_llk_call <- parLapply(cl, as.list(1:n_call), function(i) { 
    cond_llk_call <- lapply(as.list(1:n_call), function(i) { 
      
      # Create a logical of the detection history of call i
      index <- as.logical(det_hist[i, ])
      
      if (USE_BEARINGS & !USE_RL) {
        # Get bearings for call i
        bearings <- bearings_rad[i, index]
        
        #################### Derive det_probs and p. #############################
      } else if (USE_RL & !USE_BEARINGS) {
        # Get received levels for call i
        rl <- received_levels[i, index]
        
        #################### Derive det_probs and p. #############################
        
        # Get noise for call i
        c <- noise_call[i, ]
        
      } else if (USE_RL & USE_BEARINGS) {
        # Get received levels and bearings for call i
        bearings <- bearings_rad[i, index]
        rl <- received_levels[i, index]
        
        #################### Derive det_probs and p. #############################
        
        # Get noise for call i
        c <- noise_call[i, ]
        
      } else {
        #################### Derive det_probs and p. #############################
      }
      
      if (USE_RL) {
        # Get the expected rl and snr for all source lvels, grid points and detectors
        E_rl <- array(data = NA, 
                      dim = c(n_sl, n_grid, n_det),
                      dimnames = list("source_level" = source_levels, 
                                      "grid_point" = 1:n_grid, 
                                      "detector" = 1:n_det))
        
        # Get the expected snr and snr for all source lvels, grid points and detectors
        E_snr <- array(data = NA, 
                       dim = c(n_sl, n_grid, n_det),
                       dimnames = list("source_level" = source_levels, 
                                       "grid_point" = 1:n_grid, 
                                       "detector" = 1:n_det))
        
        # Create array with identical noise in every row
        c_m <- matrix(c, nrow = n_grid, ncol = n_det, byrow = TRUE)
        
        # system.time({
        for (sl in source_levels) {
          sl_m <- matrix(sl, nrow = n_grid, ncol = n_det)
          # Derive the expected snr
          E_rl[as.character(sl), , ] <- sl_m - beta_r * log(distances, base = 10)
          
          # Derive the expected snr
          E_snr[as.character(sl), , ] <- sl_m - beta_r *
            log(distances, base = 10) - c_m
        }
        
        # Set E_snr to 0 if negative
        E_snr[E_snr < 0] <- 0
        # })
        
        # Clean environment
        rm(sl, sl_m)
        
        # Derive the associated detection probabilities 
        # system.time({
        if (det_function == "janoschek" | det_function == "logit" | det_function == "probit") {
          det_probs <- .gSNR(snr = E_snr, 
                             par = par_det,
                             type = det_function)
        } else if (det_function == "half-normal") {
          # Create same distance array for every source level
          distance_array <- array(data = NA,
                                  dim = c(n_sl, n_grid, n_det),
                                  dimnames = list("source_level" = source_levels, 
                                                  "grid_point" = 1:n_grid, 
                                                  "detector" = 1:n_det))
          for (j in 1:n_sl) {distance_array[j, , ] <- distances}
          det_probs <- .gHN(distances = distance_array, 
                            par = par_det)
        }
        # })
        
        # Derive probabilities that a call was detected at least 
        # min_no_detections times
        p. <- apply(det_probs, 1, function(probs) {
          .detected(probs = probs, min_no_detections = min_no_detections)
        })
      } else {
        # Derive the associated detection probabilities 
        det_probs <- .gHN(distances = distances, 
                          par = par_det)
        p. <- .detected(probs = det_probs, min_no_detections = min_no_detections)
      }
      # Run some test whether the data is correct
      if (any(det_probs < 0)) {
        print("At least one element in det_probs is negative!")
      }
      if (any(p. < 0)) {
        print("At least one element in p. is negative!")
      }
      
      ############## Derive the different conditional likelihood parts ###########
      
      if (!USE_RL) {n_sl <- 1}
      #### Part 1: detection histories and density 
      ## Density
      part_density <- matrix(log(D), nrow = n_grid, ncol = n_sl, byrow = FALSE)
      
      ## Detection histories
      if (USE_RL) {
        part_det_hist <- apply(det_probs, 1, function(x) {
          # Create matrix of probabilities with correct detection/non-detection
          p <- 1 - x
          p[, index] <- x[, index]
          
          # Take the logarithm
          log_p <- log(p)

          # Return vector of the sum of log_p for every grid point
          return(Rfast::colsums(t(log_p)))
        })
      } else {
        part_det_hist <- apply(t(det_probs), 2, function(x) { #### SWITCH THIS TO COLUMNWISE
          # Create matrix of probabilities with correct detection/non-detection
          p <- 1 - x
          p[index] <- x[index]
          
          # Take the logarithm
          log_p <- log(p)
          
          # Return vector of the sum of log_p for every grid point
          return(sum(log_p))
        })
      }
      part_1 <- part_density + part_det_hist
      # Replace -Inf with -6000, since exp(-6000) = 0
      part_1[part_1 == -Inf] <- -6000
      
      #### Part 2: bearings
      if (USE_BEARINGS) {
        
        obs_minus_exp <- circular::circular(matrix(bearings, nrow = n_grid, ncol = length(bearings), byrow = TRUE) - 
                                              grid_bearings[, index], template = "geographics")
        
        log_p <- circular::dvonmises(x = obs_minus_exp, 
                                     mu = circular::circular(0, template = "geographics"),
                                     kappa = kappa,
                                     log = TRUE)
        # Return the sum of the log probabilities
        part_bearings <- colSums(t(log_p))
       
      } else {part_bearings <- 0}
      
      part_2 <- matrix(part_bearings, nrow = n_grid, ncol = n_sl, byrow = FALSE)
      # Replace -Inf with -6000, since exp(-6000) = 0
      part_2[part_2 == -Inf] <- -6000
      
      ## Part 3: received levels and source level !THIS IS THE SLOW PART!
      # system.time({
      if (USE_RL) {
        # part_received_levels <- 0
        
        part_received_levels <- t(apply(E_rl[, , index], c(1, 2), function(x) {
          # Subtract expected levels from received levels
          # obs_minus_exp <- rl - x
          
          # # dtruncnorm() does not allow for the use of logarithms, but very fast
          # p <- truncnorm::dtruncnorm(x = rl,
          #                            a = c[index],
          #                            mean = x,
          #                            sd = sd_r)
          # p_log <- log(p)
          
          # Assume normal on received levels instead of truncated normal
          p_log <- dnorm(x = rl, mean = x, sd = sd_r, log = TRUE)
          
          return(sum(p_log))
        }))
        
        part_source_level <- truncnorm::dtruncnorm(x = matrix(source_levels, nrow = n_grid, 
                                                              ncol = n_sl, byrow = TRUE), 
                                                   a = lower_s, b = upper_s, mean = mu_s, 
                                                   sd = sd_s)
        part_source_level <- log(part_source_level)
        
        # Add the two subpart to create part 3
        part_3 <- part_received_levels + part_source_level
      } else {
        part_3 <- matrix(0, nrow = n_grid, ncol = n_sl)
      }
      # Replace -Inf with -6000, since exp(-6000) = 0
      part_3[part_3 == -Inf] <- -6000
      # })
      
      ## Part 4: the effective sampled area given c_i, a|c_i
      if (USE_RL) {
        # Calculate log(p.) + log(D) + log(f(s))
        temp <- log(p.) + part_density + part_source_level
        
        # logSumExp() first over the source levels
        temp <- apply(t(temp), 2, matrixStats::logSumExp)
      } else {
        # Calculate log(p.) + log(D) 
        temp <- log(p.) + part_density
      }
      
      # logSumExp() over the grid points
      part_4 <- log(A_s) + log(A_x) + matrixStats::logSumExp(temp)
      # Replace -Inf with -6000, since exp(-6000) = 0
      part_4[part_4 == -Inf] <- -6000
      
      rm(temp)
      
      #### Add all the parts together
      total_per_grid_sl <- part_1 + part_2 + part_3 - part_4
      
      # Check for -Inf values due to logarithms of 0
      # replace values probabilities smaller than log(error) with log(error).
      # This might cuase problem of (almost) all p's are smaller than log(error)
      # Or use -6000 as smallest on log scale
      total_per_grid_sl[is.infinite(total_per_grid_sl) | total_per_grid_sl < -6000] <- -6000
      
      # logSumExp() over the grid points (OLD: and add log of area A)
      total_per_sl <- apply(total_per_grid_sl, 2, matrixStats::logSumExp)
      
      if (USE_RL) {
        if (smooth == "gam") {
          # Fit a smooth 
          smooth_fit <- mgcv::gam(total_per_sl ~ s(source_levels)) # gam is more robust but slower
          pred <- mgcv::predict.gam(smooth_fit,
                                    newdata = data.frame(source_levels = source_levels_all))
        } else if (smooth == "loess") {
          smooth_fit <- loess(total_per_sl ~ source_levels) # loess has edge effects
          pred <- predict(smooth_fit,
                          newdata = data.frame(source_levels = source_levels_all))
        } else if (smooth == "none") {pred <- total_per_sl}
        
        # logSumExp() over the source levels and add the logs of the areas
        total <- matrixStats::logSumExp(pred) + log(A_s) + log(A_x)
      } else {total <- total_per_sl}
      
      return(total)
    })
    
    # stopCluster(cl)
    
    # Sum all the cond llk's for every call to get the total cond llk
    llk_cond <- sum(unlist(cond_llk_call))
    
    
    ########################### Create complete llk ############################
    
    llk_full <- llk_n + llk_cond
  }
  
  ##############################################################################
  ########################## Non-parallel version ##############################
  ##############################################################################
  # system.time({
  if (!USE_RL) {
    ################ Derive the estimate for rate parameter ####################
    ## Derive probabilities of being detected enough times
    # Derive the detection probabilities based on distances
    det_probs <- .gHN(distances = distances, par = par_det)
    
    # Derive probabilities that a call was detected at least 
    # min_no_detections times for all source levels
    p. <- .detected(probs = det_probs, min_no_detections = min_no_detections)
    
    # Multiply p. and D, then sum and multiply by A to get the rate parameter
    # or ESA in case of constant density (I think?)
    rate <- A_x * sum(p. * D)
    
    if (!CONSTANT_DENSITY) {
      part_full <- dpois(n_call, lambda = rate, log = TRUE)

      if (is.infinite(part_full)) {part_full <- -6000}
    } else {part_full <- 0}
    
    ################ Derive the conditional likelihood part ####################
    #### Part 1: detection histories and density 
    ## Detection histories
    part_det_hist <- apply(det_hist, 1, function(w) {
      # Create matrix of probabilities with correct detection/non-detection
      index <- as.logical(w)
      p <- det_probs
      p[, !index] <- (1 - det_probs)[, !index]
      
      # Take the logarithm
      log_p <- log(p)
      
      # Return vector of the sum of log_p for every grid point
      return(colSums(t(log_p)))
    })
    ## Density
    part_density <- log(matrix(D, nrow = n_grid, ncol = n_call, byrow = FALSE))
    
    ### Part 2: bearings
    if (USE_BEARINGS) {
      part_bearings <- apply(bearings_rad, 1, function(x) {
        # Subtract bearings for grid from observed bearings to centre on zero
        index <- !is.na(x)
        obs_minus_exp <- circular(matrix(x, nrow = n_grid, ncol = n_det, byrow = TRUE) - 
                                    grid_bearings, template = "geographics")
        obs_minus_exp <- obs_minus_exp[, index] ## I THINK THIS IS NOT CORRECT
        log_p <- dvonmises(x = obs_minus_exp, 
                           mu = circular::circular(0, template = "geographics"),
                           kappa = kappa,
                           log = TRUE)
        
        # Return the sum of the log probabilities
        return(colSums(t(log_p)))
      })
    } else {part_bearings <- 0}
    
    # Add all parts together
    part_conditional <- part_det_hist + part_density + part_bearings
    # Some over the grid points using the logSumExp approximations and add log(A) - log(rate)
    part_conditional <- apply(part_conditional, 2, logSumExp) - log(rate) + log(A_x)
    part_conditional[is.infinite(part_conditional)] <- -6000
    # Sum over all the calls
    part_conditional <- sum(part_conditional)
    
    ####
    llk_full <- part_full + part_conditional
  }
  # })
  
  if (TRACE) {cat("The log likelihood is: ", llk_full, "\n", sep = "")}
  
  return(llk_full)
}
