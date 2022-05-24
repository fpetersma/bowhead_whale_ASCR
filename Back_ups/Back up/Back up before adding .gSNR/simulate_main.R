################################################################################
####                Simulate data for performance testing                   ####
################################################################################

## Load required libraries
library("readr")
library("matrixStats")
library("truncnorm")

source("Scripts/Function_v2/hidden_functions.R")
source("Scripts/Function_v2/plotDensity.R")
source("Scripts/Function_v2/simulateData.R")
source("Scripts/Function_v2/bwASCR.R")
source("Scripts/Function_v2/llkLSEParallelSmooth.R")

## Define some constants
seed <- 333
sample_size <- 30
min_no_detections <- 2
max_depth <- 150

## Load a fine grid of the study area
density_grid <- read.csv("Data/Masks/[PT] mask_buffer=30km_mesh=coarse_A=8774550m2_deltaX=3145_deltaY=2790_grid=421.csv")
cov_density <- density_grid
colnames(cov_density)[c(1,2)] <- c("y", "x")
cov_density$depth[cov_density$depth == 0] <- 0.01
cov_density$depth[cov_density$depth > max_depth] <- max_depth
cov_density[["depth2"]] <-  cov_density$depth ^ 2
cov_density[["logdepth"]] <- log(cov_density$depth)
cov_density[["dist_to_coast2"]] <- cov_density$dist_to_coast ^ 2

# Set the area parameter 
A <- 8.774550 # km^2

## Load the detector location data
DASAR <- as.data.frame(read_tsv("Data/DASARs.txt"))
DASAR$long <- -DASAR$long # Turn longitude from westing to easting (more common)
detectors <-  DASAR[DASAR$year == 2010 & DASAR$site == 5, ] 
detectors <- detectors[, c("utmx", "utmy")]
colnames(detectors) <- c("x", "y")

# Define the density function
f_density <- D ~ 1
# f_density <- D ~ dist_to_coast + dist_to_coast2
# f_density <- D ~ depth + depth2

# set det function!
# det_function <- "half-normal"
det_function <- "jano"

## Set all parameters on the real scale. Only density is on the log scale
# Parameters for density function
par_dens <- c("(Intercept)" = -1.5)
# par_dens <- c("(Intercept)" = -2, "dist_to_coast" = 1, "dist_to_coast2" = -1.5)
# par_dens <- c("(Intercept)" = -3, "depth" = 0.5, "depth2" = -2)

if (det_function == "jano") {
  # Detection function parameters
  par_det <- c(U = 0.8, # should be in (0, 1]
               B = 2, # should be in (0, Inf)
               Q = 4) # should be in (1, Inf)
} else {
  par_det <- c(g0 = 0.7,
               sigma = 15000)
}

# Distribution of mean noise level per call
par_noise <- c(mu = 40,
               sd = 10,
               lower = 0)

# Source level distribution
par_sl <- c(mu = 100,
            sd = 5,
            lower = 70,
            upper = 130)

# Received level parameters
par_rl <- c(beta = 15,
            sigma = 1)

# Bearing distribution
par_bear <- c(kappa = 10)

par <- list(par_dens = par_dens,
            par_det = par_det,
            par_noise = par_noise,
            par_sl = par_sl,
            par_rl = par_rl,
            par_bear = par_bear)

# Set seed
set.seed(seed)

################ Run the simulation using simulateData.R #######################
dat_sim <- simulateData(par = par, f_density = f_density, cov_density = cov_density,
                        A = A, min_no_detections = min_no_detections,
                        detectors = detectors, det_function = det_function)

################# Fit ASCR model to the simulated data #########################

sample_size <- nrow(dat_sim$det_hist) # 30

index <- sample(1:nrow(dat_sim$det_hist), sample_size)

det_hist <- dat_sim$det_hist[index, ]
bearings_hist <- dat_sim$bearings[index, ]
noise_call <- dat_sim$noise_call[index, ]
noise_random <- dat_sim$noise_random[1:50, ]
received_levels_hist <- dat_sim$received_levels[index, ]

DASAR <- as.data.frame(read_tsv("Data/DASARs.txt"))

mask <- read.csv("Data/Masks/[PT] mask_buffer=30km_mesh=coarse_A=8774550m2_deltaX=3145_deltaY=2790_grid=421.csv")

########################### Create mask ########################################

colnames(mask)[1:2] <- c("y", "x")
cov_density <- mask
cov_density$depth[cov_density$depth == 0] <- 0.01
cov_density$depth[cov_density$depth > max_depth] <- max_depth
cov_density[["depth2"]] <-  cov_density$depth ^ 2
cov_density[["logdepth"]] <- log(cov_density$depth)
cov_density[["dist_to_coast2"]] <- cov_density$dist_to_coast ^ 2

######################### Create detectors #####################################

DASAR$long <- -DASAR$long # Turn longitude from westing to easting (more common)
detectors <-  DASAR[DASAR$year == 2010 & DASAR$site == 5, ]
detectors <- detectors[, c("utmx", "utmy")]
colnames(detectors) <- c("x", "y")

# Density formula specification
f_density <- D ~ 1
# f_density <- D ~ dist_to_coast + dist_to_coast2
# f_density <- D ~ depth + depth2

## Define some constants
min_no_detections <- 2

# Set the area parameter
A <- 8.774550 # km^2

det_function <- "jano"
# det_function <- "half-normal"

dat <- list(det_hist = det_hist,
            detectors = detectors,
            cov_density = cov_density,
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
if(det_function == "jano") {
  par_det_start <- c(U = log(0.5 / (1 - 0.5)), # should be in (0, 1]
                     B = log(1), # should be in (0, Inf)
                     Q = log(5 - 1)) # should be in (1, Inf)
} else {
  par_det_start <- c(g0 = log(0.7 / (1 - 0.7)),
                     sigma = log(15000))
}

# Start values for density function
par_dens_start <- c("(Intercept)" = -0)
# par_dens_start <- c("(Intercept)" = -1.5, "dist_to_coast" = 1, "dist_to_coast2" = -1.5)
# par_dens_start <- c("(Intercept)" = -1, "depth" = -1, "depth2" = 3)

# Start values for received level
par_rl_start <- c(beta_r = log(5),
                  sd_r = log(3))

# Start values for source level
par_sl_start <- c(mu_s = log(80), # identity
                  sd_s = log(10)) # log

# Start values for bearing
par_bear_start <- c(kappa = log(20)) # log

par_start <- list(par_det = par_det_start,
                  par_bear = par_bear_start,
                  par_rl = par_rl_start,
                  par_sl = par_sl_start,
                  par_dens = par_dens_start)

######################## Run bwASCR() on subset data ##########################
method = "L-BFGS-B"
maxit = 100
TRACE = TRUE
LSE = TRUE

system.time({
  bw <- bwASCR(dat = dat, par = par_start, method = method, maxit = maxit,
               TRACE = TRACE, LSE = LSE)
})

