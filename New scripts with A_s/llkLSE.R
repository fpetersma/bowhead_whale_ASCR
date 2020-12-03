# for test runs
par <- par 
dat <- dat
rm(list = setdiff(ls(), c("dat", "par")))



# The 'fast' loglikelihood
.llkLSE <- function(par, dat) {
  
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

  if (TRACE) {print(par)}
  
  if (USE_BEARINGS) {
    bearings_rad <- dat$bearings_rad
    grid_bearings <- dat$grid_bearings
  }
  if (USE_RL) {
    source_levels <- dat$source_levels
    received_levels <- dat$received_levels
    snr <- dat$SNR
    noise <- dat$noise
  }
  
  ################## Extract some important constants ##########################
  n <- nrow(det_hist) # number of calls
  K <- nrow(detectors) # number of detectors
  M <- nrow(cov_density) # number of grid points
  
  vars <- colnames(design) # parameters for density function
  
  # Extract parameters and convert to the real scale (add errors to ensure correct domain)
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
  
  ################## Derive probabilites for all combinations ##################
  
  # Derive the density
  D <- .density_GAM(x = design, gam_fit = gam_fit, par = par_dens)
  
  # Derive expected snr for all source levels, grid points, calls and detectors
  E_snr <- array(data = NA, 
                 dim = c(length(source_levels), n, M, K),
                 dimnames = list("source_level" = source_levels, 
                                 "call" = 1:n, 
                                 "grid_point" = 1:M, 
                                 "detector" = 1:K))
  det_probs <- E_snr
  system.time({
    for (sl in source_levels) { # Check the runtime of this hideous monster
      # print(sl)
      for (i in 1:n) {
        c <- matrix(noise[i, ], nrow = M, ncol = K, byrow = TRUE)
        # Derive the expected snr
        E_snr[as.character(sl), i, , ] <- sl - beta_r * log(distances, base = 10) - c
        
        # # Derive the associated detection probabilities (doing this here is slooooow)
        # det_probs[as.character(sl), i, , ] <- .p(snr = E_snr[as.character(sl), i, , ],
        #                                          par = par_det,
        #                                          det_function = det_function)
        # 
        # # Derive the associated probability of being detected at least the specified no of times
        # # so slloooooooooooooooooow
        # p.[as.character(sl), i, ] <- .detected(det_probs[as.character(sl), i, , ],
        #                                        min_no_detections = min_no_detections)
      }
    }
  })
  if (any(!is.finite(E_snr))) {stop()} # Check for non finite values
  rm(sl, c, i)
  
  # Derive detection probabilities (doing this this way is really fast)
  system.time({
    det_probs <- .p(snr = E_snr, par = par_det, det_function = det_function)
  })
  
  # Derive probabilities that a call was detected at least min_no_detections times
  p. <- array(data = NA,
              dim = c(length(source_levels), n, M),
              dimnames = list("source_level" = source_levels,
                              "call" = 1:n,
                              "grid_point" = 1:M))
  ### Using foreach() similar to a for loop #######################
  library(parallel)
  library(foreach)
  library(doParallel)
  
  system.time({
  no_cores <- detectCores() - 1 # get number of cores - 1
  registerDoParallel(no_cores) # initiate parallelisation
  # Use foreach() to run parallel operations on the data for every source level 
  p._list <- foreach(sl = 70:130) %dopar% {
    temp <- matrix(nrow = n, ncol = M)
    for (i in 1:n) {
        # p.[as.character(sl), i, ] <- .detected(probs = det_probs[as.character(sl), i, , ],
        #                                          min_no_detections = min_no_detections)
        temp[i, ] <- .detected(probs = det_probs[as.character(sl), i, , ],
                                  min_no_detections = min_no_detections)
    }
    return(temp)
    # t(sapply(1:n, function(i) {
    #   .detected(probs = a[i, , ], min_no_detections = min_no_detections)
    # }))
  }
  # clean up the cluster
  stopImplicitCluster()
  # Add list elements of p._list to the correct array locations in p.
  for (i in seq_along(source_levels)) {p.[i, , ] <- p._list[[i]]} 
  })
  
  ### Creating a list of one dimension and than parallel lapply() ##############
  # # Faster than foreach()
  # library(parallel)
  # library(snow)
  # 
  # system.time({
  # det_probs_list <- lapply(seq(dim(det_probs)[1]), function(x) det_probs[x, , ,])
  # 
  # no_cores <- detectCores() - 1
  # cl <- makeCluster(no_cores)
  # clusterExport(cl, list("min_no_detections", ".detected", "n", "M"))
  # p._list <- parLapply(cl, det_probs_list, function(a) {
  #   # temp <- matrix(nrow = n, ncol = M)
  #   # for (i in 1:n) {
  #   #   temp[i, ] <- .detected(probs = a[i, , ], min_no_detections = min_no_detections)
  #   # }
  #   # return(temp)
  #   
  #   temp <- t(sapply(1:n, function(i) {
  #     .detected(probs = a[i, , ], min_no_detections = min_no_detections)
  #   }))
  #   return(temp)
  # })
  # stopCluster(cl)
  # # Add list elements of p._list to the correct array locations in p.
  # for (i in seq_along(source_levels)) {p.[i, , ] <- p._list[[i]]} 
  # })
  ##############################################################################

  ########## This is some old code that works correctly ########################
  # 
  # system.time({
  # for (sl in source_levels) { # Check the runtime of this hideous monster
  #   print(sl)
  #   for (i in 1:n) {
  #     # for(m in 1:M) {
  #       # p.[as.character(sl), i, ] <-  rep(1, M)#.detected(det_probs = det_probs[as.character(sl), i, , ],
  # 
  #       # microbenchmark({                               #        min_no_detections = min_no_detections)
  #       p.[as.character(sl), i, ] <- .detected(probs = det_probs[as.character(sl), i, , ],
  #                                              min_no_detections = min_no_detections)
  #       # })
  #     # }
  #   }
  # }
  # })
  # 
  # system.time({
  # a <- apply(det_probs, c(1, 2), function(probs) {
  #   .detected(probs = probs, min_no_detections = min_no_detections)
  # })
  # })
  ##############################################################################
  
  # Run some test whether the data is correct
  if (any(D < 0)) {
    print("At least one element in D is negative!")
  }
  if (any(det_probs < 0)) {
    print("At least one element in det_probs is negative!")
  }
  if (any(p. < 0)) {
    print("At least one element in p. is negative!")
  }
  
  ## PART 1: A *
  #         sum of density times probability of detected more than twice /
  #         modified bessel function of degree 0 of kappa to the power of n times K
  # f_x and p. are already derived, just derive f_s
  
  f_s <- dnorm(source_levels, mu_s, sigma_s) # normal
  # f_s <- dtruncnorm(source_levels, a = 70, b = 130, mean = 100, sd = 10) # truncated norm
  part_1a <- p. * array(f_s, dim = dim(p.)) * aperm(array(D, dim = rev(dim(p.))))
  
  part_1 <- -A * sum(D * p.) # * n * K * log(besselI(kappa, 0))
  
  ## Part 2
  if (USE_BEARINGS) {
    # Limit kappa to 1000000 (resulting in near perfect detection)
    if (kappa > log(100000)) {
      kappa_real <- 100000
    } else {
      kappa_real <- exp(kappa)
    }
    
    part_b <- log(2 * pi) + log(besselI(kappa_real, 0, TRUE)) + kappa_real
    
    part_2_bearings_log <- apply(as.matrix(bearings_rad), 1, function(y) {
      index <- !is.na(y)
      
      ########## Calculate the probabilities for the bearings ####################
      y_matrix <- matrix(y[index], nrow = M, ncol = sum(index), byrow = TRUE)
      differences <- y_matrix - as.matrix(grid_bearings[, index])
      
      part_a <- kappa_real * cos(differences)
      # sum for an individual for every grid point
      part_bearings_im <- rowSums(part_a - part_b)
      
      return(part_bearings_im)
    }) 
  } else {
    part_2_bearings_log <- 0
  }
  part_2_det_log <- apply(as.matrix(det_hist), 1, function(w) {
    ######### Calculate the probabilities for detections #####################
    index <- as.logical(w)
    p <- matrix(NA, nrow = M, ncol = K)
    p[, index] <- det_probs[, index]
    p[, !index] <- (1 - det_probs)[, !index]
    
    if (any(is.na(p))) {
      stop("Error in calculation of part_2_2, NA values detected!")
    }
    
    p[p < 1e-16] <- 1e-16
    log_p <- log(p)
    
    # Take the sum of the log of all probabilities for every individual at every gridpoint
    p_im <- rowSums(log_p)
    
    return(p_im)
  })
  # Part 2: the signal strength data
  if (USE_SS) {
    part_2_s_log <- list()
    sigma_s_real <- exp(sigma_s)
    b0_real <- exp(b0)
    b1_real <- exp(b1)
    r_expected <- b0_real - b1_real * log(distances, base = 10) # Sometimes r_expected becomes 0, should I set these to 0 or is this fine?
    signal_strength <- as.matrix(signal_strength)
    noise <- as.matrix(noise)
    
    for (i in 1:nrow(signal_strength)) {
      index <- !is.na(signal_strength[i, ])
      r_i <- matrix(signal_strength[i, index], nrow = M, ncol = sum(index), byrow = TRUE) # signal strength for call i at all detectors
      c_i <- matrix(noise[i, index], nrow = M, ncol = sum(index), byrow = TRUE) # noise for call i at all detectors
      r_i_expected <- r_expected[, index]
      
      a <- dnorm((r_i - r_i_expected) / sigma_s_real, log = TRUE)
      b <- log(1 - pnorm((c_i - r_i_expected) / sigma_s_real) + 1e-16)
      
      prob_r_ik <- -matrix(log(sigma_s_real), nrow = M, ncol = sum(index)) + a - b
      prob_r_im <- rowSums(prob_r_ik)
      part_2_s_log[[i]] <- prob_r_im
    }
    part_2_s_log <- as.matrix(data.frame(part_2_s_log))
    colnames(part_2_s_log) <- 1:n
  } else {part_2_s_log <- matrix(0, nrow = M, ncol = n)}
  # part_2_ss_log <- apply(as.matrix(det_hist), 1, function(w) {})
  
  # Calculate a matrix of the sum of all log probailities for individual i and gridpoint m
  part_2im_log <- t(log(D) + part_2_det_log +  part_2_bearings_log + part_2_s_log)
  if (nrow(part_2im_log) != n | ncol(part_2im_log) != M) {
    stop("Error in calculation of part_2_ik, dimensions are incorrect!")
  }
  
  # Take the log-summmation over M-exp using logSumExp()
  part_2i_log <- apply(part_2im_log, 1, logSumExp)
  if (any(is.infinite(part_2i_log))) {
    warning("Error in calculation of part_2_i, log(-x) not possible. Try different starting values.")
  }
  
  # Return the sum of the logs of the likelihoods of every call
  part_2 <- sum(part_2i_log)
  
  # cat(paste0("part1:", part_1, "\n"))
  # cat(paste0("part2:", part_2, "\n\n"))
  
  ## Part 3: modelling the signal-to-noise process
  if (USE_SNR) {
    sigma_s_real <- exp(sigma_s)
    b0_real <- exp(b0)
    b1_real <- exp(b1)
    
    SNR_expected <- b0_real - b1_real * log(distances, base = 10)
    SNR_expected[SNR_expected <= 0] <- 1e-16 # Set really small number to avoid negative SNR. Maybe use a log link function?
    
    part_SNR_log <- t(apply(as.matrix(SNR), 1, function(y) {
      index <- !is.na(y)
      
      ########## Calculate the probabilities for the bearings ##################
      y_matrix <- matrix(y[index], nrow = M, ncol = sum(index), byrow = TRUE)
      differences <- log(y_matrix) - log(as.matrix(SNR_expected[, index]))
      
      # part_a <- 1 / (y_matrix * sigma_s_real * sqrt(2 * pi)) * exp(- (differences ^ 2) / (2 * sigma_s_real ^ 2))
      part_a_log <- 0 - log(y_matrix * sigma_s_real * sqrt(2 * pi)) - (differences ^ 2) / (2 * sigma_s_real ^ 2)
      
      # product for an individual for every grid point
      # part_SNR_im <- rowProds(part_a) # Not using log for stability
      part_SNR_im <- rowSums(part_a_log) # Using log which results in sums for stability
      
      return(part_SNR_im)
    }))
    # Make sure the dimensions are correct (rows are the individuals, columns the gridpoints)
    if (nrow(part_SNR_log) != n | ncol(part_SNR_log) != M) {
      stop("Error in calculation of part_SNR, dimensions are incorrect!")
    }
    # Now take the log-sum-exp for every individuals over all gridpoints
    part_3i_log <- apply(part_SNR_log, 1, logSumExp)
    
    part_3 <- sum(part_3i_log)
  } else {
    part_3 <- 0
  }
  
  total <- part_1 + part_2 + part_3
  cat(total)
  cat("\n")
  
  if (total == -Inf) {total <- -1e100}
  if (total == Inf) {cat("A positive return? That's wierd..\n")}
  
  return(total)
}