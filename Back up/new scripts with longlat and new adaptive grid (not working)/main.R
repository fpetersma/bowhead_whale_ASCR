################################################################################
####                              "main.R"                                  ####
####  --------------------------------------------------------------------  ####
####  This script contains  all the required functionality for the fitting. ####
####                                                                        ####
#### ---------------------------------------------------------------------- ####
#### By:        Felix Thomas Petersma                                       ####
#### Purpose:   Developed as part of the MSc Statistics dissertation at     ####
####            the University of St Andrews [2018/2019]                    ####
#### License:   MIT                                                         ####
################################################################################
library(readr)
library(microbenchmark)

# run the script that sources all functions
source("Scripts/Function_v2/hidden_functions.R")
source("Scripts/Function_v2/llkLSE_parallel.R")
source("Scripts/Function_v2/bwASCR.R")

########################### Define constants ###################################

max_depth <- 150 # Declare the max depth for the depth covariate to improve stability
seed <- 121
sample_size <- 100
min_no_detections <- 2

############################ Load real data ####################################

detections <- read.csv("Data/Revised data with different cut-offs [22-11-19]/[PT] detection_history_S5Y10_0.02.csv")
bearings <- read.csv("Data/Revised data with different cut-offs [22-11-19]/[PT] bearings_history_S5Y10_0.02.csv")
weights <- read.csv("Data/Revised data with different cut-offs [22-11-19]/[PT] weights_history_S5Y10_0.02.csv")
received_levels <- read.csv("Data/Revised data with different cut-offs [22-11-19]/[PT] sound_levels_history_S5Y10_0.02.csv")
snr <- read.csv("Data/Revised data with different cut-offs [22-11-19]/[PT] SNR_history_S5Y10_0.02.csv")

DASAR <- as.data.frame(read_tsv("Data/DASARs.txt"))

mask <- read.csv("Data/Masks/[PT] mask_buffer=75km_mesh=coarse_A=26862864m2_deltaX=5158_deltaY=5208_grid=406.csv")

########################### Create mask ########################################

colnames(mask)[1:2] <- c("y", "x")
cov_density <- mask
cov_density$depth[cov_density$depth == 0] <- 0.01
# cov_density$depth[cov_density$depth > max_depth] <- 150
cov_density[["depth2"]] <-  cov_density$depth ^ 2
cov_density[["logdepth"]] <- log(cov_density$depth)
cov_density[["dist_to_coast2"]] <- cov_density$dist_to_coast ^ 2

######################### Create detectors #####################################

DASAR$long <- -DASAR$long # Turn longitude from westing to easting (more common)
detectors <-  DASAR[DASAR$year == 2010 & DASAR$site == 5, ] 
detectors <- detectors[, c("utmx", "utmy")]
colnames(detectors) <- c("x", "y")

######################### Create usuable data ##################################
# Remove all detections with weight below 0.9
weights[is.na(weights) | weights < 0.9] <- 0
weights[weights != 0] <- 1
bearings[weights == 0] <- NA
detections[weights == 0] <- 0

# Filter out linkages with too few detections THIS SHOULD BE DONE OUTSIDE OF THIS SCRIPT
too_few_detections <- rowSums(!is.na(bearings)) < min_no_detections
detections <- detections[!too_few_detections, ]
bearings <- bearings[!too_few_detections, ]
received_levels <- received_levels[!too_few_detections, ]
snr <- snr[!too_few_detections, ]

# Set seed for reproducibility
# set.seed(seed)

# Detection history
index <- sample(1:nrow(detections), sample_size)
det_hist <- detections[row.names(detections)[index], ]
bearings_hist <- bearings[row.names(bearings)[index], ]
snr_hist <- snr[row.names(bearings)[index], ]
received_levels_hist <- received_levels[row.names(bearings)[index], ]
noise_call <- received_levels_hist - snr_hist
# If noise contains NA, replace those by the mean of the noise. 
noise_call <- t(apply(noise_call, 1, function(c) {
  k <- is.na(c)
  c[k] <- rep(mean(c, na.rm = TRUE), sum(k)) + rnorm(sum(k), 0, 2) # added noise 
                                                                   # to noise
  return(c)
}))

# Create fake random sample of noise based on noise_call
noise_random <- noise_call + rnorm(length(noise_call), 4, 2)
noise_random <- noise_random[sample(1:nrow(noise_random), 30), ]

sample_percentage <- nrow(det_hist) / nrow(detections)

# Density formula specification
f_density <- D ~ 1 #s(dist_to_coast, k = 3, fx = TRUE) # + depth
# f_density <- D ~ dist_to_coast + dist_to_coast2
# f_density <- D ~ s(depth, k = 3, fx = TRUE)
# f_density <- D ~ dist_to_coast + dist_to_coast2 + depth + depth2


# Area per grid point (in square kilometers)
A <- 26.862864 # make sure this is correct for the mask!

det_function <- "jano"
det_function <- "half-normal"

dat <- list(det_hist = det_hist, 
            detectors = detectors,
            cov_density = cov_density, 
            # snr = snr_hist, 
            received_levels = received_levels_hist,
            noise_call = noise_call,
            noise_random = noise_random,
            source_levels = 70:130,
            A = A, 
            f_density = f_density,
            det_function = det_function, 
            bearings = bearings_hist, 
            min_no_detections = min_no_detections)

# Start values for detection function
par_det_start <- c(U = 0, # should be in (0, 1]
                   B = 0, # should be in (0, Inf)
                   Q = 0) # should be in (1, Inf)
par_det_start <- c(g0 = 0.5,
                   sigma = log(30000))

# Start values for density function
par_dens_start <- c("(Intercept)" = -2.5)

# Start values for received level
par_rl_start <- c(beta_r = log(15),
                  sd_r = log(1))

# Start values for source level
par_sl_start <- c(mu_s = log(100),
                  sd_s = log(5))

# Start values for bearing
par_bear_start <- c(kappa = log(25))

par_start <- list(par_det = par_det_start, 
                  par_dens = par_dens_start, 
                  par_rl = par_rl_start,
                  par_sl = par_sl_start,
                  par_bear = par_bear_start)

######################## Run bw_scr() on subset data ##########################
method = "L-BFGS"
maxit = 100
TRACE = TRUE
LSE = TRUE

system.time({
  bw <- bwASCR(dat = dat, par = par_start, method = method, maxit = maxit, 
                TRACE = TRACE, LSE = LSE)
})
