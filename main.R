# ============================================================================ #
#                                 "main.R"                                     #
# ---------------------------------------------------------------------------- #
#     This script contains  all the required functionality for the fitting.    #
#                                                                              #
# ---------------------------------------------------------------------------- #
# By:        Felix Thomas Petersma                                             #
# Purpose:   Developed as part of the MSc Statistics dissertation at           #
#            the University of St Andrews [2018/2019]                          #
# License:   MIT                                                               #
# ============================================================================ #

# Load required libraries
library("readr")
library("matrixStats")
library("truncnorm")

# run the script that sources all functions
source("Scripts/bowhead whales/hidden_functions.R")
source("Scripts/bowhead whales/plotDensity.R")
source("Scripts/bowhead whales/simulateData.R")
source("Scripts/bowhead whales/bwASCR.R")
source("Scripts/bowhead whales/llkParallelSmooth.R")

# source("Scripts/Bowhead Whales/bwASCRNoSL.R")
# source("Scripts/Bowhead Whales/llkLSEParallelSmoothNoSL.R")

########################### Define constants ###################################

# max_depth <- 150 # Declare the max depth for the depth covariate to improve stability
seed <- 31082010
sample_size <- 500
min_no_detections <- 2
SINGLE_SL <- TRUE
WITH_NOISE <- FALSE
trunc_level <- 96
BEAR_MIXTURE <- TRUE

############################ Load real data ####################################

detections <- read.csv("Data/detections_31-08-2010_all.csv")
bearings <- read.csv("Data/bearings_31-08-2010_all.csv")
received_levels <- read.csv("Data/received_levels_31-08-2010_all.csv")
noise_call <- read.csv("Data/noise_31-08-2010_all.csv") # 99% quantile of noise is 95.6

DASAR <- as.data.frame(read_tsv("Data/DASARs.txt"))

mesh_file <- "Data/grid_adaptive_levels=1_maxD2C=100k_maxD2A=100k_area=21785.36_n=1163.csv"

# Create fake random sample of noise based on noise_call
noise_random <- noise_call + rnorm(length(noise_call), 2, 1)
noise_random <- noise_random[sample(1:nrow(noise_random), 100), ]

########################### Create mask ########################################

mesh <- read.csv(mesh_file)
cov_density <- mesh
cov_density$depth <- abs(cov_density$depth) # make depth positive 
# cov_density$depth[cov_density$depth == 0] <- 0.01
cov_density[["depth2"]] <-  cov_density$depth ^ 2
cov_density[["logdepth"]] <- log(cov_density$depth)
cov_density[["logdepth2"]] <- cov_density$logdepth ^ 2
cov_density[["distance_to_coast2"]] <- cov_density$distance_to_coast ^ 2

######################### Create detectors #####################################

DASAR$long <- -abs(DASAR$long) # Turn longitude from westing to easting (more common)
detectors <- DASAR[DASAR$year == 2010 & DASAR$site == 5, c("long", "lat")] 

######################### Create usable data ##################################
# # Remove all detections with weight below 0.9
# weights[is.na(weights) | weights < 0.9] <- 0
# weights[weights != 0] <- 1
# bearings[weights == 0] <- NA
# detections[weights == 0] <- 0

# Filter out linkages with too few detections THIS SHOULD BE DONE OUTSIDE OF THIS SCRIPT
# too_few_detections <- rowSums(!is.na(bearings)) < min_no_detections
# detections <- detections[!too_few_detections, ]
# bearings <- bearings[!too_few_detections, ]
# received_levels <- received_levels[!too_few_detections, ]
# snr <- snr[!too_few_detections, ]


#===========# Run next section to filter on received levels #==================#

if (WITH_NOISE) {
  sufficient_snr <- received_levels - noise_call >= trunc_level # get indices 
  sufficient_snr[is.na(sufficient_snr)] <- FALSE
  
  detections <- sufficient_snr * 1
  enough_dets <- Rfast::rowsums(detections) > 1
  
  # only keep calls with enough detections, after truncation at 15dB snr
  sufficient_snr <- sufficient_snr[enough_dets, ]
  detections <- detections[enough_dets, ]
  
  bearings <- as.matrix(bearings)[enough_dets, ]
  bearings[!sufficient_snr] <- NA
  
  noise_call <- as.matrix(noise_call)[enough_dets, ]
  
  received_levels <- as.matrix(received_levels)[enough_dets, ]
  received_levels[!sufficient_snr] <- NA
} else if (!WITH_NOISE) {
  sufficient_rl <- received_levels  >= trunc_level # get indices 
  sufficient_rl[is.na(sufficient_rl)] <- FALSE
  
  detections <- sufficient_rl * 1
  enough_dets <- Rfast::rowsums(detections) > 1
  
  # only keep calls with enough detections, after truncation at 15dB snr
  sufficient_rl <- sufficient_rl[enough_dets, ]
  detections <- detections[enough_dets, ]
  
  bearings <- as.matrix(bearings)[enough_dets, ]
  bearings[!sufficient_rl] <- NA
  
  noise_call <- as.matrix(noise_call)[enough_dets, ]
  
  received_levels <- as.matrix(received_levels)[enough_dets, ]
  received_levels[!sufficient_rl] <- NA
}

#==============================================================================#

# Set seed for reproducibility
set.seed(seed)

# Detection history
# index <- sample(1:nrow(detections), sample_size)
det_hist <- detections#[index, ]
bearings_hist <- bearings#[index, ]
received_levels_hist <- received_levels#[index, ]
noise_call_hist <- noise_call#[index, ]

# get sample percentage
sample_percentage <- nrow(det_hist) / nrow(detections)

# Density formula specification
f_density <- D ~ 1 #s(dist_to_coast, k = 3, fx = TRUE) # + depth
# f_density <- D ~ s(distance_to_coast, k = 5, fx = TRUE)
# f_density <- D ~ distance_to_coast + distance_to_coast2
# f_density <- D ~ s(logdepth, k = 3, fx = TRUE) + s(distance_to_coast, k = 3, fx = TRUE)
# f_density <- D ~ dist_to_coast + dist_to_coast2 + depth + depth2


# Set the source level integration increments and the grid areas
A_s <- 3
A_x <- subset(mesh, select = c(area, long, lat))

det_function <- "simple"

dat <- list(det_hist = as.matrix(detections),
            detectors = as.matrix(detectors),
            cov_density = as.data.frame(cov_density),
            received_levels = as.matrix(received_levels_hist),
            noise_call = as.matrix(noise_call),
            noise_random = as.matrix(noise_random),
            source_levels = seq(from = 150 + (A_s / 2), to = 200, by = A_s),
            A_x = A_x,
            A_s = A_s,
            f_density = f_density,
            det_function = det_function,
            bearings = as.matrix(bearings_hist),
            min_no_detections = min_no_detections,
            SINGLE_SL = SINGLE_SL,
            WITH_NOISE = WITH_NOISE,
            BEAR_MIXTURE = BEAR_MIXTURE,
            trunc_level = trunc_level)


# Start values for detection function
if(det_function == "janoschek") {
  par_det_start <- c(U = log(0.8 / (1 - 0.8)), # should be in (0, 1]
                     B = log(1.5), # should be in (0, Inf)
                     Q = log(2.8 - 1)) # should be in (1, Inf)
} else if (det_function == "logit" | det_function == "probit") {
  # Detection function parameters
  par_det_start <- c(U = log(0.7 / (1 - 0.7)), # should be in (0, 1]
                     B = log(10), # should be in (0, Inf)
                     Q = log(10)) # should be in (0, Inf)
} else if (det_function == "simple") { ## STILL NEEDS IMPLEMENTATION AND TRUNCATION
  par_det_start <- c(g0 = log(0.6 / (1 - 0.6)))
} else {
  par_det_start <- c(g0 = log(0.8 / (1 - 0.8)),
                     sigma = log(15000))
}

# Start values for density function
par_dens_start <- c("(Intercept)" = -3)
# par_dens_start <- c("(Intercept)" = -3, "distance_to_coast" = 2, "distance_to_coast2" = -2.5)
# par_dens_start <- c("(Intercept)" = -3, "depth" = 0.5, "depth2" = -2)

# Start values for received level
par_rl_start <- c(beta_r = log(15),
                  sd_r = log(4))

# Start values for source level
par_sl_start <- c(mu_s = log(154), # identity
                  sd_s = log(8)) # log

# Start values for bearing
if (BEAR_MIXTURE) {
  par_bear_start <- c(kappa_low = log(1), #log
                      kappa_high = log(10),
                      mix_par = log(0.1 / (1 - 0.1))) # logit
} else {
  par_bear_start <- c(kappa = log(5)) #log
} 

par_start <- list(par_det = par_det_start,
                  par_bear = par_bear_start,
                  par_rl = par_rl_start,
                  par_sl = par_sl_start,
                  par_dens = par_dens_start)



######################## Run bw_scr() on subset data ##########################
method = "L-BFGS-B"
maxit = 100
TRACE = 1
LSE = TRUE


system.time({
  bw <- bwASCR(dat = dat, par = par_start, method = method, maxit = maxit,
               TRACE = TRACE, LSE = LSE)
})


