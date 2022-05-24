# ============================================================================ #
# Simulation script to test the performance of the four variations on the      #
# method. 20 data set will be simulated with noise and variable source         #
# levels, and then all four models will be fitted to these data sets.          #
# ============================================================================ #

# Load libraries 
library("readr")
library("matrixStats")
library("truncnorm")
library("tidyverse")

# source the essential scripts
source("Scripts/Bowhead Whales/hidden_functions.R")
source("Scripts/Bowhead Whales/plotDensity.R")
source("Scripts/Bowhead Whales/simulateData.R")
source("Scripts/Bowhead Whales/bwASCR.R")
source("Scripts/Bowhead Whales/llkParallelSmooth.R")

# ==============================================================================
# Simulate the 20 data sets                                                    
# ------------------------------------------------------------------------ Start
seeds <- seq(from = 090596, length.out = 20)

mesh_file <- "Data/grid_adaptive_levels=3_bounds=10k-60k_maxD2C=100k_maxD2A=200k_area=44145.6_n=180+139+146=465.csv"

simulated_data <- lapply(seeds, function(seed) {
  # declare constants
  min_no_detections <- 2
  trunc_level <- 10
  SINGLE_SL <- TRUE
  WITH_NOISE <- TRUE
  
  # load a fine grid of the study area
  mesh <- read.csv(mesh_file)
  cov_density <- mesh
  cov_density$depth <- abs(cov_density$depth) # make depth positive 
  cov_density$depth[cov_density$depth == 0] <- 0.01
  cov_density[["depth2"]] <-  cov_density$depth ^ 2
  cov_density[["logdepth"]] <- log(cov_density$depth)
  cov_density[["logdepth2"]] <- cov_density$logdepth ^ 2
  cov_density[["distance_to_coast2"]] <- cov_density$distance_to_coast ^ 2
  
  # load the detector location data
  DASAR <- as.data.frame(read_tsv("Data/DASARs.txt"))
  DASAR$long <- -abs(DASAR$long) # Turn longitude from westing to easting (more common)
  detectors <- DASAR[DASAR$year == 2010 & DASAR$site == 5, c("long", "lat")] 
  
  # define the density function
  f_density <- D ~ distance_to_coast + distance_to_coast2

  # set detection function
  det_function <- "simple"
  if (det_function == "simple") par_det <- c(g0 = 0.7)
  
  # set all parameters on the real scale. Only density is on the log scale
  # parameters for density function
  par_dens <- c("(Intercept)" = -2, 
                "distance_to_coast" = 3, 
                "distance_to_coast2" = -3)
  
  # distribution of mean noise level per call
  par_noise <- c(mu = 65,
                 sd = 2,
                 lower = 0,
                 upper = Inf)
  
  # source level distribution 
  par_sl <- c(mu = 140,
              sd = 5,
              lower = 0,
              upper = Inf)
  
  # received level parameters
  par_rl <- c(beta = 15,
              sigma = 3)
  
  # bearing distribution
  par_bear <- c(kappa = 100)
  
  par <- list(par_dens = par_dens,
              par_det = par_det,
              par_noise = par_noise,
              par_sl = par_sl,
              par_rl = par_rl,
              par_bear = par_bear)
  
  # set seed
  set.seed(seed)
  
  # run the simulation using simulateData()
  dat_sim <- simulateData(par = par, 
                          f_density = f_density, 
                          cov_density = cov_density,
                          min_no_detections = min_no_detections, 
                          SINGLE_SL = SINGLE_SL, 
                          WITH_NOISE = WITH_NOISE,
                          detectors = detectors, 
                          det_function = det_function,
                          trunc_level = trunc_level)
  
  return(dat_sim)
})
# -------------------------------------------------------------------------  End

# ==============================================================================
# Create function that truncates data                                      Begin 
# ------------------------------------------------------------------------------
truncateData <- function(d, trunc_level, min_detections = 2) {
  # extract received levels, truncate and find calls with sufficient detections
  rl <- d$received_levels
  rl[is.na(rl)] <- 0
  rl[rl < trunc_level] <- NA
  sufficient_det <- rowSums(!is.na(rl)) >= min_detections
  
  # filter all but noise_random in d
  truncated_data <- lapply(d[-4], function(x) {
    return(x[sufficient_det, ])
  })
  
  # add noise_random
  truncated_data[["noise_random"]] <- d$noise_random 
  
  return(truncated_data)
}
# -------------------------------------------------------------------------- End


# ==============================================================================
# Single source level without noise                                        Begin
# ------------------------------------------------------------------------------
system.time({
fits_singleSL_rl <- lapply(1:20, function(i) {
  
  # extract simulated data
  dat_sim <- simulated_data[[i]]
  
  # set constants
  SINGLE_SL <- TRUE
  WITH_NOISE <- FALSE
  trunc_level <- 75 # which is the 100% quantile = ~75 (the 99.99% quantile = ~72)
  
  # filter data due to truncation
  dat_sim <- truncateData(d = dat_sim, 
                          trunc_level = trunc_level, 
                          min_detections = 2)
  print(nrow(dat_sim$det_hist))
  # # extract data from dat_sim
  # det_hist <- dat_sim$det_hist
  # bearings_hist <- dat_sim$bearings
  # noise_call <- dat_sim$noise_call
  # noise_random <- dat_sim$noise_random[1:100, ] # get a noise sample of size 100
  # received_levels_hist <- dat_sim$received_levels
  # 
  # DASAR <- as.data.frame(read_tsv("Data/DASARs.txt"))
  # 
  # mesh <- read.csv(mesh_file)
  # 
  # ########################### Create mask ########################################
  # cov_density <- mesh
  # cov_density$depth <- abs(cov_density$depth) # make depth positive 
  # cov_density$depth[cov_density$depth == 0] <- 0.01
  # cov_density[["depth2"]] <-  cov_density$depth ^ 2
  # cov_density[["logdepth"]] <- log(cov_density$depth)
  # cov_density[["logdepth2"]] <- cov_density$logdepth ^ 2
  # cov_density[["distance_to_coast2"]] <- cov_density$distance_to_coast ^ 2
  # 
  # ######################### Create detectors #####################################
  # 
  # DASAR$long <- -DASAR$long # Turn longitude from westing to easting (more common)
  # detectors <-  DASAR[DASAR$year == 2010 & DASAR$site == 5, c("long", "lat")]
  # 
  # # Density formula specification
  # f_density <- D ~ distance_to_coast + distance_to_coast2
  # 
  # # Define some constants
  # min_no_detections <- 2
  # 
  # # Set the source level integration increments and the grid areas
  # A_s <- 3
  # A_x <- subset(mesh, select = c(area, long, lat))
  # 
  # # Set detection function
  # det_function <- "simple"
  # 
  # dat <- list(det_hist = det_hist,
  #             detectors = detectors,
  #             cov_density = cov_density,
  #             received_levels = received_levels_hist,
  #             noise_call = noise_call,
  #             noise_random = noise_random,
  #             source_levels = seq(from = 120 + (A_s / 2), to = 160, by = A_s),
  #             A_x = A_x,
  #             A_s = A_s,
  #             f_density = f_density,
  #             det_function = det_function,
  #             bearings = bearings_hist,
  #             min_no_detections = min_no_detections,
  #             SINGLE_SL = SINGLE_SL,
  #             WITH_NOISE = WITH_NOISE,
  #             trunc_level = trunc_level)
  # 
  # if (det_function == "simple") {
  #   par_det_start <- c(g0 = log(0.6 / (1 - 0.6)))
  # }
  # 
  # # Start values for density function
  # par_dens_start <- c("(Intercept)" = -2, 
  #                     "distance_to_coast" = 3, 
  #                     "distance_to_coast2" = -3)
  # 
  # # Start values for received level
  # par_rl_start <- c(beta_r = log(15),
  #                   sd_r = log(3))
  # 
  # # Start values for source level
  # par_sl_start <- c(mu_s = log(140), # identity
  #                   sd_s = log(5)) # log
  # 
  # # Start values for bearing
  # par_bear_start <- c(kappa = log(20)) # log
  # 
  # par_start <- list(par_det = par_det_start,
  #                   par_bear = par_bear_start,
  #                   par_rl = par_rl_start,
  #                   par_sl = par_sl_start,
  #                   par_dens = par_dens_start)
  # 
  # ######################## Run bw_scr() on subset data ##########################
  # method = "L-BFGS-B"
  # maxit = 100
  # TRACE = 0
  # LSE = TRUE
  # 
  # system.time({
  #   bw <- bwASCR(dat = dat, par = par_start, method = method, maxit = maxit,
  #                TRACE = TRACE, LSE = LSE)
  # })
  # 
  # save(list = c("bw"), 
  #      file = paste0("save singleSL rl index=", i, ".RData"))
  # 
  # return(bw)
})
})
# -------------------------------------------------------------------------- End


# ==============================================================================
# Single source level with noise                                           Begin
# ------------------------------------------------------------------------------
fits_singleSL_snr <- lapply(1:20, function(index) {
  # extract simulated data
  dat_sim <- simulated_data[[i]]
  
  # extract data from dat_sim
  det_hist <- dat_sim$det_hist
  bearings_hist <- dat_sim$bearings
  noise_call <- dat_sim$noise_call
  noise_random <- dat_sim$noise_random[1:100, ] # get a noise sample of size 100
  received_levels_hist <- dat_sim$received_levels
  
  DASAR <- as.data.frame(read_tsv("Data/DASARs.txt"))
  
  mesh <- read.csv(mesh_file)
  
  ########################### Create mask ########################################
  cov_density <- mesh
  cov_density$depth <- abs(cov_density$depth) # make depth positive 
  cov_density$depth[cov_density$depth == 0] <- 0.01
  cov_density[["depth2"]] <-  cov_density$depth ^ 2
  cov_density[["logdepth"]] <- log(cov_density$depth)
  cov_density[["logdepth2"]] <- cov_density$logdepth ^ 2
  cov_density[["distance_to_coast2"]] <- cov_density$distance_to_coast ^ 2
  
  ######################### Create detectors #####################################
  
  DASAR$long <- -DASAR$long # Turn longitude from westing to easting (more common)
  detectors <-  DASAR[DASAR$year == 2010 & DASAR$site == 5, c("long", "lat")]
  
  # Density formula specification
  f_density <- D ~ distance_to_coast + distance_to_coast2
  
  # Define some constants
  min_no_detections <- 2
  
  # Set the source level integration increments and the grid areas
  A_s <- 3
  A_x <- subset(mesh, select = c(area, long, lat))
  
  # Set detection function
  det_function <- "simple"
  
  dat <- list(det_hist = det_hist,
              detectors = detectors,
              cov_density = cov_density,
              received_levels = received_levels_hist,
              noise_call = noise_call,
              noise_random = noise_random,
              source_levels = seq(from = 120 + (A_s / 2), to = 180, by = A_s),
              A_x = A_x,
              A_s = A_s,
              f_density = f_density,
              det_function = det_function,
              bearings = bearings_hist,
              min_no_detections = min_no_detections,
              SINGLE_SL = SINGLE_SL,
              WITH_NOISE = WITH_NOISE,
              trunc_level = trunc_level)
  
  if (det_function == "simple") {
    par_det_start <- c(g0 = log(0.7 / (1 - 0.7)))
  }
  
  # Start values for density function
  par_dens_start <- c("(Intercept)" = -2.5, 
                      "distance_to_coast" = 3, 
                      "distance_to_coast2" = -3)
  
  # Start values for received level
  par_rl_start <- c(beta_r = log(15),
                    sd_r = log(3))
  
  # Start values for source level
  par_sl_start <- c(mu_s = log(140), # identity
                    sd_s = log(5)) # log
  
  # Start values for bearing
  par_bear_start <- c(kappa = log(10)) # log
  
  par_start <- list(par_det = par_det_start,
                    par_bear = par_bear_start,
                    par_rl = par_rl_start,
                    par_sl = par_sl_start,
                    par_dens = par_dens_start)
  
  ######################## Run bw_scr() on subset data ##########################
  method = "L-BFGS-B"
  maxit = 100
  TRACE = 0
  LSE = TRUE
  
  system.time({
    bw <- bwASCR(dat = dat, par = par_start, method = method, maxit = maxit,
                 TRACE = TRACE, LSE = LSE)
  })
  
  save(list = c("bw"), 
       file = paste0("temp_saves/save singleSL snr index=", i, ".RData"))
  
  return(bw)
})
# -------------------------------------------------------------------------- End


# ==============================================================================
# Variable source level without noise 
# ------------------------------------------------------------------------ Start
fits_variableSL_rl
# -------------------------------------------------------------------------- End

# ==============================================================================
# Variable source level with noise 
# ------------------------------------------------------------------------ Start
fits_variableSL_snr
# -------------------------------------------------------------------------- End
