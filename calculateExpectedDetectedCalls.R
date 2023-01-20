## =============================================================================
## calculatedExpectedDetectedCalls.R
##
## This script derives the expected number of call detection histories that 
## involved only 1 sensor for model 33 (and other scenarios, for fun)
## =============================================================================

## =============================================================================
## 1. LOAD DATA AND LIBRARIES
## -----------------------------------------------------------------------------
load("Real data output/fits_1_35_nlminb_n=443.RData")

fit <- fits[[33]]

pars <- fit$real_par
dat <- fit$dat

# Remove the last row which contains the calibration for no detection per call
min_det_matrix <- dat$min_det_matrix[-7, ]

## =============================================================================
## 2. SET PARAMETERS
## -----------------------------------------------------------------------------
# Loop over grid cell and store densities corrected for detection prob and f(s)
results <- rep(NA, nrow(fit$D))

for (grid_no in 1:nrow(fit$D)) {
  print(grid_no)
  # Get distances to sensors
  distances <- dat$X[grid_no, ]

  source_level_support <- 100:220
  
  # Loop over source level and store results
  corrected_det_probs_per_sl <- rep(NA, length(source_level_support))
  for (sl_index in seq_along(source_level_support)) {
    sl <- source_level_support[sl_index]
    
    # Derive expected received levels at all detectors
    E_rl <- sl - pars["log_beta_r"] * log10(distances)
    
    # Derive positive and negative detection probabilities
    p1 <- matrix(pars["logit_g0"] * (1 - pnorm((dat$trunc_level - E_rl) / pars["log_sd_r"])),
                  nrow = 6, ncol = 6, byrow = TRUE)
    
    # # Version without uncertainty on the received level
    # p1 <- matrix(ifelse(dat$trunc_level > E_rl, pars["logit_g0"], 0), 
    #              nrow = 6, ncol = 6, byrow = TRUE)
    
    p0 <- 1 - p1
    
    # Combine those in single matrix for all combinations that give 1 positive detection
    det_probs <- p0
    det_probs[min_det_matrix == 1] <- p1[min_det_matrix == 1]
    
    # Derive total probability through row products and a sum
    p.1 <- sum(Rfast::rowprods(det_probs))
    
    # p.1(x, s) * f(s)
    corrected_det_probs_per_sl[sl_index] <- p.1 * dnorm(sl, 
                                              mean = pars["log_mu_s"], 
                                              sd = pars["log_sd_s"])
  }
  total_p.1xs_fs <- sum(corrected_det_probs_per_sl)
  
  # Store number of singletons per grid cell in results
  results[grid_no] <- fit$D$N[grid_no] * total_p.1xs_fs
}

## Total number of singletons = ~587.36
total_singletons <- sum(results)


## =============================================================================
## Total number of calls with > 1 detections 
min_2_det_calls <- nrow(dat$W)
true_frequencies <- table(Rfast::rowsums(dat$W))

## Ratio expected singletons vs observed singletons (= )
## Load the real data
detections <- read.csv("Data/detections_31-08-2010_all.csv")
received_levels <- read.csv("Data/received_levels_31-08-2010_all.csv")

sufficient_rl <- received_levels  >= dat$trunc_level # get indices 
sufficient_rl[is.na(sufficient_rl)] <- FALSE

detections <- sufficient_rl * 1
singletons <- Rfast::rowsums(detections) == 1

# only keep calls with enough detections, after truncation at 96db RL
sufficient_rl <- sufficient_rl[singletons, ]
detections <- detections[singletons, ]

## 723 observed singletons from the min 2 detection data 
##    + 770 from the original singleton data!! = 1493 singletons



## =============================================================================
# Expected number of zerotons, just for shits ==================================
min_det_matrix <- dat$min_det_matrix[7, ] # + 1 # the + 1 is to calculate all positive detections
# Loop over grid cell and store densities corrected for detection prob and f(s)
results <- rep(NA, nrow(fit$D))

for (grid_no in 1:nrow(fit$D)) {
  print(grid_no)
  # Get distances to detetors
  distances <- dat$X[grid_no, ]
  
  # Loop over source level and store results
  corrected_det_probs_per_sl <- rep(NA, 121)
  for (sl_index in seq_along(100:220)) {
    sl <- (100:220)[sl_index]
    
    # Derive expected received levels at all detectors
    E_rl <- sl - pars["log_beta_r"] * log10(distances)
    
    # Derive positive and negative detection probabilities
    p1 <- matrix(pars["logit_g0"] * (1 - pnorm((dat$trunc_level - E_rl) / pars["log_sd_r"])),
                 nrow = 1, ncol = 6, byrow = TRUE)
    p0 <- 1 - p1
    
    # Combine those in single matrix for all combinations that no positive detections
    det_probs <- p0
    det_probs[min_det_matrix == 1] <- p1[min_det_matrix == 1]
    
    # Derive total probability through row products and a sum
    p.1 <- sum(Rfast::rowprods(det_probs))
    
    # p.1(x, s) * f(s)
    corrected_det_probs_per_sl[sl_index] <- p.1 * dnorm(sl, 
                                                        mean = pars["log_mu_s"], 
                                                        sd = pars["log_sd_s"])
  }
  total_p.1xs_fs <- sum(corrected_det_probs_per_sl)
  
  # Store number of calls with 1 detection per grid cell in results
  results[grid_no] <- fit$D$N[grid_no] * total_p.1xs_fs
}

sum(results) # ~4739, seems to add up

## =============================================================================
## Try to find estimated for twice detected calls

ff <- as.matrix(expand.grid(rep(list(0:1), 6)))
min_det_matrix <- ff[rowsums(ff) == 2, ]

# Loop over grid cell and store densities corrected for detection prob and f(s)
results <- rep(NA, nrow(fit$D))

for (grid_no in 1:nrow(fit$D)) {
  print(grid_no)
  # Get distances to sensors
  distances <- dat$X[grid_no, ]
  
  source_level_support <- 100:220
  
  # Loop over source level and store results
  corrected_det_probs_per_sl <- rep(NA, length(source_level_support))
  for (sl_index in seq_along(source_level_support)) {
    sl <- source_level_support[sl_index]
    
    # Derive expected received levels at all detectors
    E_rl <- sl - pars["log_beta_r"] * log10(distances)
    
    # Derive positive and negative detection probabilities
    p1 <- matrix(pars["logit_g0"] * (1 - pnorm((dat$trunc_level - E_rl) / pars["log_sd_r"])),
                 nrow = nrow(min_det_matrix), ncol = 6, byrow = TRUE)
    
    # # Version without uncertainty on the received level
    # p1 <- matrix(ifelse(dat$trunc_level > E_rl, pars["logit_g0"], 0), 
    #              nrow = 6, ncol = 6, byrow = TRUE)
    
    p0 <- 1 - p1
    
    # Combine those in single matrix for all combinations that give 1 positive detection
    det_probs <- p0
    det_probs[min_det_matrix == 1] <- p1[min_det_matrix == 1]
    
    # Derive total probability through row products and a sum
    p.1 <- sum(Rfast::rowprods(det_probs))
    
    # p.1(x, s) * f(s)
    corrected_det_probs_per_sl[sl_index] <- p.1 * dnorm(sl, 
                                                        mean = pars["log_mu_s"], 
                                                        sd = pars["log_sd_s"])
  }
  total_p.1xs_fs <- sum(corrected_det_probs_per_sl)
  
  # Store number of singletons per grid cell in results
  results[grid_no] <- fit$D$N[grid_no] * total_p.1xs_fs
}

## Total number of twice detected calls  = ~256.77
total_doubletons <- sum(results)


## =============================================================================
## Try to find estimated for trice detected calls

ff <- as.matrix(expand.grid(rep(list(0:1), 6)))
min_det_matrix <- ff[rowsums(ff) == 3, ]

# Loop over grid cell and store densities corrected for detection prob and f(s)
results <- rep(NA, nrow(fit$D))

for (grid_no in 1:nrow(fit$D)) {
  print(grid_no)
  # Get distances to sensors
  distances <- dat$X[grid_no, ]
  
  source_level_support <- 100:220
  
  # Loop over source level and store results
  corrected_det_probs_per_sl <- rep(NA, length(source_level_support))
  for (sl_index in seq_along(source_level_support)) {
    sl <- source_level_support[sl_index]
    
    # Derive expected received levels at all detectors
    E_rl <- sl - pars["log_beta_r"] * log10(distances)
    
    # Derive positive and negative detection probabilities
    p1 <- matrix(pars["logit_g0"] * (1 - pnorm((dat$trunc_level - E_rl) / pars["log_sd_r"])),
                 nrow = nrow(min_det_matrix), ncol = 6, byrow = TRUE)
    
    # # Version without uncertainty on the received level
    # p1 <- matrix(ifelse(dat$trunc_level > E_rl, pars["logit_g0"], 0), 
    #              nrow = 6, ncol = 6, byrow = TRUE)
    
    p0 <- 1 - p1
    
    # Combine those in single matrix for all combinations that give 1 positive detection
    det_probs <- p0
    det_probs[min_det_matrix == 1] <- p1[min_det_matrix == 1]
    
    # Derive total probability through row products and a sum
    p.1 <- sum(Rfast::rowprods(det_probs))
    
    # p.1(x, s) * f(s)
    corrected_det_probs_per_sl[sl_index] <- p.1 * dnorm(sl, 
                                                        mean = pars["log_mu_s"], 
                                                        sd = pars["log_sd_s"])
  }
  total_p.1xs_fs <- sum(corrected_det_probs_per_sl)
  
  # Store number of singletons per grid cell in results
  results[grid_no] <- fit$D$N[grid_no] * total_p.1xs_fs
}

## Total number of trice detected calls  = ~129.29
total_tripletons <- sum(results)


## =============================================================================
## Try to find estimated for quadriple detected calls

ff <- as.matrix(expand.grid(rep(list(0:1), 6)))
min_det_matrix <- ff[rowsums(ff) == 4, ]

# Loop over grid cell and store densities corrected for detection prob and f(s)
results <- rep(NA, nrow(fit$D))

for (grid_no in 1:nrow(fit$D)) {
  print(grid_no)
  # Get distances to sensors
  distances <- dat$X[grid_no, ]
  
  source_level_support <- 100:220
  
  # Loop over source level and store results
  corrected_det_probs_per_sl <- rep(NA, length(source_level_support))
  for (sl_index in seq_along(source_level_support)) {
    sl <- source_level_support[sl_index]
    
    # Derive expected received levels at all detectors
    E_rl <- sl - pars["log_beta_r"] * log10(distances)
    
    # Derive positive and negative detection probabilities
    p1 <- matrix(pars["logit_g0"] * (1 - pnorm((dat$trunc_level - E_rl) / pars["log_sd_r"])),
                 nrow = nrow(min_det_matrix), ncol = 6, byrow = TRUE)
    
    # # Version without uncertainty on the received level
    # p1 <- matrix(ifelse(dat$trunc_level > E_rl, pars["logit_g0"], 0), 
    #              nrow = 6, ncol = 6, byrow = TRUE)
    
    p0 <- 1 - p1
    
    # Combine those in single matrix for all combinations that give 1 positive detection
    det_probs <- p0
    det_probs[min_det_matrix == 1] <- p1[min_det_matrix == 1]
    
    # Derive total probability through row products and a sum
    p.1 <- sum(Rfast::rowprods(det_probs))
    
    # p.1(x, s) * f(s)
    corrected_det_probs_per_sl[sl_index] <- p.1 * dnorm(sl, 
                                                        mean = pars["log_mu_s"], 
                                                        sd = pars["log_sd_s"])
  }
  total_p.1xs_fs <- sum(corrected_det_probs_per_sl)
  
  # Store number of singletons per grid cell in results
  results[grid_no] <- fit$D$N[grid_no] * total_p.1xs_fs
}

## Total number of quadriple detected calls  = ~55.41
total_quadriple <- sum(results)

