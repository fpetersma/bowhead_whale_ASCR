################################################################################
####                Simulate data for performance testing                   ####
################################################################################

## Load required libraries
library("readr")
library("matrixStats")
library("truncnorm")

SINGLE_SL <- TRUE

source("Scripts/Bowhead Whales/hidden_functions.R")
source("Scripts/Bowhead Whales/plotDensity.R")
source("Scripts/Bowhead Whales/simulateData.R")

if (SINGLE_SL) {
  source("Scripts/Bowhead Whales/bwASCRNoSL.R")
  source("Scripts/Bowhead Whales/llkLSEParallelSmoothNoSL.R")
} else {
  source("Scripts/Bowhead Whales/bwASCR.R")
  source("Scripts/Bowhead Whales/llkLSEParallelSmooth.R")
}

mesh_file <- "Data/JABES paper/grid_adaptive_levels=3_bounds=10k-60k_maxD2C=100k_maxD2A=200k_area=44145.6_n=180+139+146=465.csv"

## Define some constants
seed <- 12326
# sample_size <- 30
min_no_detections <- 2
# max_depth <- 150

## Load a fine grid of the study area
mesh <- read.csv(mesh_file)
cov_density <- mesh
cov_density$depth <- abs(cov_density$depth) # make depth positive 
cov_density$depth[cov_density$depth == 0] <- 0.01
cov_density[["depth2"]] <-  cov_density$depth ^ 2
cov_density[["logdepth"]] <- log(cov_density$depth)
cov_density[["logdepth2"]] <- cov_density$logdepth ^ 2
cov_density[["distance_to_coast2"]] <- cov_density$distance_to_coast ^ 2

# # Set the area parameter DIFFERENT NOW AS AREA IS VARIABLE
# A <- 8.774550 # km^2
# cov_density$area <- 56.21183

## Load the detector location data
DASAR <- as.data.frame(read_tsv("Data/DASARs.txt"))
DASAR$long <- -abs(DASAR$long) # Turn longitude from westing to easting (more common)
detectors <- DASAR[DASAR$year == 2010 & DASAR$site == 5, c("long", "lat")] 

# Define the density function
f_density <- D ~ 1
# f_density <- D ~ distance_to_coast + distance_to_coast2
# f_density <- D ~ depth + depth2

# set det function!
det_function<- "probit"
# det_function <- "half-normal"
# det_function <- "logit"
# det_function <- "janoschek"

## Set all parameters on the real scale. Only density is on the log scale
# Parameters for density function
par_dens <- c("(Intercept)" = -2)
# par_dens <- c("(Intercept)" = -2.5, "distance_to_coast" = 3, "distance_to_coast2" = -3)
# par_dens <- c("(Intercept)" = -3, "depth" = 0.5, "depth2" = -2)
# par_dens <- c("(Intercept)" = -2.2, "depth" = -1.8, "depth2" = 1.8)

if (det_function == "janoschek") {
  # Detection function parameters
  par_det <- c(U = 0.8, # should be in (0, 1]
               B = 1.5, # should be in (0, Inf)
               Q = 2.8) # should be in (1, Inf)
} else if (det_function == "logit" | det_function == "probit") {
  # Detection function parameters
  par_det <- c(U = 0.7, # should be in (0, 1]
               B = 20, # should be in (0, Inf)
               Q = 5) # should be in (0, Inf)
} else {
  par_det <- c(g0 = 0.7,
               sigma = 15000)
}

# Distribution of mean noise level per call
par_noise <- c(mu = 65,
               sd = 8,
               lower = 0,
               upper = Inf)

# Source level distribution 
par_sl <- c(mu = 140,
            sd = 5,
            lower = 0,
            upper = Inf)

# Received level parameters
par_rl <- c(beta = 15,
            sigma = 5)

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
                        min_no_detections = min_no_detections, 
                        SINGLE_SL = SINGLE_SL,
                        detectors = detectors, det_function = det_function)

################# Fit ASCR model to the simulated data #########################
# sample_size <- nrow(dat_sim$det_hist) # 30
index <- #1:nrow(dat_sim$det_hist) #sample(1:nrow(dat_sim$det_hist), sample_size)

det_hist <- dat_sim$det_hist#[index, ]
bearings_hist <- dat_sim$bearings#[index, ]
noise_call <- dat_sim$noise_call#[index, ]
noise_random <- dat_sim$noise_random[1:200, ]
received_levels_hist <- dat_sim$received_levels#[index, ]

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
f_density <- D ~ 1 #s(dist_to_coast, k = 3, fx = TRUE) # + depth
# f_density <- D ~ distance_to_coast + distance_to_coast2

# Define some constants
min_no_detections <- 2

# Set the source level integration increments and the grid areas
A_s <- 3
A_x <- subset(mesh, select = c(area, long, lat))
# A_x$area <- 56.21183

det_function <- "probit"
# det_function <- "janoschek"
# det_function <- "logit"
# det_function <- "half-normal"

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
            min_no_detections = min_no_detections)


# Start values for detection function
if(det_function == "janoschek") {
  par_det_start <- c(U = log(0.8 / (1 - 0.8)), # should be in (0, 1]
                     B = log(1.5), # should be in (0, Inf)
                     Q = log(2.8 - 1)) # should be in (1, Inf)
} else if (det_function == "logit" | det_function == "probit") {
  # Detection function parameters
  par_det_start <- c(U = log(0.7 / (1 - 0.7)), # should be in (0, 1]
                     B = log(20), # should be in (0, Inf)
                     Q = log(5)) # should be in (0, Inf)
} else {
  par_det_start <- c(g0 = log(0.8 / (1 - 0.8)),
                     sigma = log(15000))
}

# Start values for density function
par_dens_start <- c("(Intercept)" = -2)
# par_dens_start <- c("(Intercept)" = -2.5, "distance_to_coast" = 3, "distance_to_coast2" = -3)
# par_dens_start <- c("(Intercept)" = -3, "depth" = 0.5, "depth2" = -2)

# Start values for received level
par_rl_start <- c(beta_r = log(15),
                  sd_r = log(5))

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
TRACE = 6
LSE = TRUE

system.time({
  bw <- bwASCR(dat = dat, par = par_start, method = method, maxit = maxit,
               TRACE = TRACE, LSE = LSE)
})

