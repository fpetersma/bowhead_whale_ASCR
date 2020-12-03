# for test runs
par <- par
dat <- dat
rm(list = setdiff(ls(), c("dat", "par")))

# The 'fast' loglikelihood (not really, but faster than otherwise)
.llkLSE_parallel <- function(par, dat) {
  
  # Description:
  #   
  
  # Inputs:
  #   param       - [list] 
  #   dat         - [list] 
  
  # Outputs:
  #   total       - [scalar] the total log likelihood value.
  
  ###################### INPUT CHECKS STILL MISSING ############################
  
  
  # Define very small error to be added to values that are numerically zero
  error <- 1e-16
  
  ######################## Extract data from dat ###############################
  det_hist <- dat$det_hist
  cov_density <- dat$cov_density
  detectors <- dat$detectors
  A <- dat$A
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
    source_levels <- dat$source_levels
    received_levels <- dat$received_levels
    snr <- dat$snr
    noise_call <- dat$noise_call
    if (!CONSTANT_DENSITY) {
      noise_random <- dat$noise_random
    }
  }
  
  ################## Extract some important constants ##########################
  n_call <- nrow(det_hist) # number of calls
  n_det <- nrow(detectors) # number of detectors
  n_grid <- nrow(cov_density) # number of grid points
  if (USE_RL) {
    n_sl <- length(source_levels) # number of source levels
    if (!CONSTANT_DENSITY) {
      n_noise <- nrow(noise_random) # number of random noise samples
    }
  }
  
  vars <- colnames(design) # parameters for density function
  
  # Extract parameters and convert to the real scale (add errors to ensure 
  # correct domain)
  if (det_function == "jano") {
    U <- exp(par["U"]) / (1 + exp(par["U"])) # logit link 
    B <- exp(par["B"]) + error # log link
    Q <- exp(par["Q"]) + 1 + error # log link + 1
    par_det <- c(U, B, Q)
  } else {stop("A PROBLEM!")}
  
  if (USE_BEARINGS) {kappa <- exp(par["kappa"])} # log link
  
  if (USE_RL) {
    beta_r <- exp(par["beta_r"]) # log link
    sigma_r <- exp(par["sigma_r"]) # log link
    
    mu_s <- par["mu_s"] # identity
    sigma_s <- exp(par["sigma_s"]) # log link
  }
  
  # Check if all density parameters have starting values
  if (!all(vars %in% names(par))) {
    stop("Not all density parameters have starting values.")
  }
  par_dens <- par[vars]
  
  ##############################################################################
  ###################### Start deriving the log likelihood #####################
  ##############################################################################
  
  ################ Derive probabilites for all combinations ####################
  
  # Derive the density
  D <- .density_GAM(x = design, gam_fit = gam_fit, par = par_dens)
  if (any(D < 0)) {
    print("At least one element in D is negative!")
  }
  
  ################ Start parallelisation over source level #####################
  
  ### Creating a list of one dimension and than parallel lapply() ##############
  # Faster than foreach()
  library(parallel)
  library(snow)

  system.time({
    if (!CONSTANT_DENSITY) {
      ##########################################################################
      ################ Derive the estimate for rate parameter ##################
      ##########################################################################
      
      ## Initiate parallel process leaving one core free
      no_cores <- detectCores() - 1
      cl <- makeCluster(no_cores)
      
      # Export required data and functions to all clusters
      hidden_functions <- c(".gJano", ".detected", ".density_GAM")
      clusterExport(cl, list = c(ls(), hidden_functions), envir = environment()) 

      # Turn noise random in a t(data frame) and then into a list to create a 
      # list with the noise for six detectors as a vector in every element.
      noise_random_list <- as.list(as.data.frame(t(noise_random)))
      
      # Start parallel process over noise_random
      rates_noise <- parLapply(cl, noise_random_list, function(c) {
      
        # Get the expected SNR for all noise samples, grid points and detectors
        E_snr <- array(data = NA,
                       dim = c(n_sl, n_grid, n_det),
                       dimnames = list("source_level" = source_levels, 
                                       "grid_point" = 1:n_grid, 
                                       "detector" = 1:n_det))
        # Create matrix with identical noise in every row
        c_m <- matrix(c, nrow = n_grid, ncol = n_det, byrow = TRUE)
        system.time({
          for (sl in source_levels) {
            sl_m <- matrix(sl, nrow = n_grid, ncol = n_det)
            # Derive the expected snr
            E_snr[as.character(sl), , ] <- sl_m - beta_r * 
              log(distances, base = 10) - c_m
          }
        })
        # Clean environment
        rm(sl_m, sl) 
        
        # Derive the associated detection probabilities 
        system.time({
          det_probs <- .gJano(snr = E_snr, 
                              par = par_det, 
                              det_function = det_function)
        })
        # Derive probabilities that a call was detected at least 
        # min_no_detections times for all source levels
        system.time({
          p. <- t(apply(det_probs, 1, function(probs) {
            .detected(probs = probs, min_no_detections = min_no_detections)
          }))
        })
        
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
        probs_sl <- dnorm(source_levels, mean = mu_s, sd = sigma_s)
        probs_sl_m <- matrix(probs_sl, nrow = n_sl, ncol = n_grid, 
                             byrow = FALSE)
        
        # Create a matrix of all products of D, p. and f(s)
        out <- p. * D_m * probs_sl_m
        
        # Sum over all columns and rows
        out <- sum(colSums(temp)) 
        
        return(out)
      }) 
    # Stop parallel process
    stopCluster(cl)
    
    # Sum the output together and multiply by A to get the rate parameter 
    rate <- A * sum(unlist(rates_noise))
    
    # Derive the probability density for the number of detections given the rate
    llk_n <- dpois(n_call, rate, log = TRUE)
    # microbenchmark(llk_n <- n_call * log(rate) - rate - log(factorial(n_call)))
    
    # Clear environment
    rm("rate", "cl", "rates_noise", "no_cores", "hidden_functions", 
       "noise_random_list")
    } else {llk_n <- 0}
  })
  ##############################################################################
  ############### Derive the conditional likelihood part #######################
  ##############################################################################
  
  ## Initiate parallel process leaving one core free
  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores)
  
  # Export required data and functions to all clusters
  hidden_functions <- c(".gJano", ".detected", ".density_GAM")
  # package_functions <- c("dvonmises")
  clusterExport(cl, list = c(ls(), hidden_functions), envir = environment()) 
  
  # # Turn noise random in a t(data frame) and then into a list to create a 
  # # list with the noise for six detectors as a vector in every element.
  # noise_call_list <- as.list(as.data.frame(t(noise_call)))
  
  # Start parallel process over number of calls 
  cond_llk_call <- parLapply(cl, as.list(1:n_call), function(i) { 
    
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
      # Get the expected rl and snr for all calls, grid points and detectors
      E_rl <- array(data = NA, 
                    dim = c(n_sl, n_grid, n_det),
                    dimnames = list("source_level" = source_levels, 
                                    "grid_point" = 1:n_grid, 
                                    "detector" = 1:n_det))
      
      # Get the expected snr and snr for all calls, grid points and detectors
      E_snr <- array(data = NA, 
                     dim = c(n_sl, n_grid, n_det),
                     dimnames = list("source_level" = source_levels, 
                                     "grid_point" = 1:n_grid, 
                                     "detector" = 1:n_det))
      
      # Create array with identical noise in every row
      c_m <- matrix(c, nrow = n_grid, ncol = n_det, byrow = TRUE)
      
      system.time({
        for (sl in source_levels) {
          sl_m <- matrix(sl, nrow = n_grid, ncol = n_det)
          # Derive the expected snr
          E_rl[as.character(sl), , ] <- sl_m - beta_r * log(distances, base = 10)
          
          # Derive the expected snr
          E_snr[as.character(sl), , ] <- sl_m - beta_r *
            log(distances, base = 10) - c_m
        }
      })
      
      # Clean environment
      rm(sl, sl_m)
      
      # Derive the associated detection probabilities 
      system.time({
        det_probs <- .gJano(snr = E_snr, 
                            par = par_det, 
                            det_function = det_function)
      })
      
      # Derive probabilities that a call was detected at least 
      # min_no_detections times
      system.time({
        p. <- apply(det_probs, 1, function(probs) {
          .detected(probs = probs, min_no_detections = min_no_detections)
        })
      })
    } else {
      # Derive the associated detection probabilities 
      system.time({
        det_probs <- .gHN(distances = distances, 
                                   par = par_det)
      })
    }
    # Run some test whether the data is correct
    if (any(det_probs < 0)) {
      print("At least one element in det_probs is negative!")
    }
    if (any(p. < 0)) {
      print("At least one element in p. is negative!")
    }
    
    ############## Derive the different conditional likelihood parts ###########
    
    #### Part 1: detection histories and density 
    ## Density
    part_density <- matrix(log(D), nrow = n_grid, ncol = n_sl, byrow = FALSE)
    
    ## Detection histories
    part_det_hist <- apply(det_probs, 1, function(x) {
      # Create matrix of probabilities with correct detection/non-detection
      p <- 1 - x
      p[, index] <- x[, index]
      
      # Take the logarithm
      log_p <- log(p)
      if(any(is.infinite(log_p))) {
        stop("In part 1: log_p contains (-)Inf values.")
      }
      
      # Return vector of the sum of log_p for every grid point
      return(colSums(t(log_p)))
    })
    
    part_1 <- part_density + part_det_hist
    
    #### Part 2: bearings
    if (USE_BEARINGS) {
      part_bearings <- apply(grid_bearings, 1, function(x) {
        # Subtract bearings for grid from observed bearings to centre on zero
        obs_minus_exp <- bearings - x[index]
        log_p <- circular::dvonmises(x = obs_minus_exp, 
                                     mu = circular::circular(0, template = "geographics"),
                                     kappa = kappa,
                                     log = TRUE)
        
        # Return the sum of the log probabilities
        return(sum(log_p))
      })
    } else {part_bearings <- 0}
    
    part_2 <- matrix(part_bearings, nrow = n_grid, ncol = n_sl, byrow = FALSE)
    
    ## Part 3: received levels and source level !THIS IS THE SLOW PART!
    system.time({
    if (USE_RL) {
      part_received_levels <- 0
      # part_received_levels <- t(apply(E_rl[, , index], c(1, 2), function(x) {
      #   # Subtract expected levels from received levels
      #   obs_minus_exp <- rl - x
      # 
      #   # dtruncnorm() does not allow for the use of logarithms, but very fast
      #   p <- truncnorm::dtruncnorm(x = rl,
      #                              a = c[index],
      #                              mean = x,
      #                              sd = sigma_r)
      #   p_log <- log(p)
      # 
      #   # # More stable version with logarithms applied where possible,
      #   # # but much slower (800 microsecs vs 18 microsecs)
      #   # pdf_part <- dnorm(x = as.numeric((rl - x) / sigma_r), log = TRUE)
      #   # cdf_part <- log(sigma_r * (1 - pnorm((c[index] - x) / sigma_r)))
      #   #
      #   # p_log <- pdf_part - cdf_part
      # 
      #   return(sum(p_log))
      # }))
      
      # # slightly less elegant but clearer maybe, roughly same runtime
      # system.time({
      # part_received_levels <- apply(E_rl[, , index], 1, function(x) {
      #   apply(x, 1, function(y) {
      #     # Subtract expected levels from received levels
      #     obs_minus_exp <- rl - y
      #     
      #     # dtruncnorm() does not allow for the use of logarithms, but very fast
      #     p <- truncnorm::dtruncnorm(x = rl,
      #                                a = c[index], 
      #                                mean = y,
      #                                sd = sigma_r)
      #     p_log <- log(p)
      #     
      #     # # More stable version with logarithms applied where possible,
      #     # # but much slower (800 microsecs vs 18 microsecs)
      #     # pdf_part <- dnorm(x = as.numeric((rl - x) / sigma_r), log = TRUE)
      #     # cdf_part <- log(sigma_r * (1 - pnorm((c[index] - x) / sigma_r)))
      #     # 
      #     # p_log <- pdf_part - cdf_part
      #     
      #     return(sum(p_log))
      #   })
      # })
      # })
      part_source_level <- dnorm(x = matrix(source_levels, nrow = n_grid, 
                                            ncol = n_sl, byrow = TRUE), 
                                 mean = mu_s,
                                 sd = sigma_s,
                                 log = TRUE)
      
      # Add the two subpart to create part 3
      part_3 <- part_received_levels + part_source_level
    } else {
      part_3 <- matrix(0, nrow = n_grid, ncol = n_sl)
    }
    })
    ## Part 4: the effective sampled area given c_i
    # Calculate log(p.) + log(D) + log(f(s))
    temp <- log(p.) + part_density + part_source_level
    
    # logSumExp() first over the source levels
    temp <- apply(t(temp), 2, matrixStats::logSumExp)

    part_4 <- A * matrixStats::logSumExp(temp)
    
    rm(temp)
    
    #### Add all the parts together
    total <- part_1 + part_2 + part_3 + part_4
    
    # Check for -Inf values due to logarithms of 0
    # replace values probabilities smaller than log(error) with log(error).
    # This might cuase problem of (almost) all p's are smaller than log(error)
    # Or use -6000 as smallest on log scale
    total[is.infinite(total) | total < -6000] <- -6000
    
    # logSumExp() over the grid points and add area A
    total <- apply(total, 2, matrixStats::logSumExp) + A
    
    # logSumExp() over the source levels
    total <- matrixStats::logSumExp(total)
    
    return(total)
  })
  
  stopCluster(cl)
  
  # Sum all the cond llk's for every call to get the total cond llk
  llk_cond <- sum(unlist(cond_llk_call))
  
  
  ########################### Create complete llk ##############################
  
  llk <- llk_n + llk_cond
  
  ## Maybe some checks on llk?
  
  return(llk)
  ##############################################################################
}
