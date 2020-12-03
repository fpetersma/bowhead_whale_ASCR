#### Exploring the likelihood surface for varying combinations of Q and B. #####

### Step 1: simulate data using simulate_main.R
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
# par_dens <- c("(Intercept)" = -2.2, "depth" = -1.8, "depth2" = 1.8)

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

################################################################################
### Step 2: loop over values and store likelihoods
#############################################################################
B_options <- seq(from =1.5, to = 5.5, by = 0.5)
Q_options <- seq(from =1.5, to = 5.5, by = 0.5)
llk_evals <- matrix(NA, nrow = 9, ncol = 9)
colnames(llk_evals) <- B_options
rownames(llk_evals) <- Q_options

llk_evals <- data.frame(B = NA, Q = NA, llk = NA)

for (B_loop in B_options) {
  for (Q_loop in Q_options) {
    sample_size <- nrow(dat_sim$det_hist) # 30
    
    index <- sample(1:nrow(dat_sim$det_hist), sample_size)
    
    det_hist <- dat_sim$det_hist[index, ]
    bearings_hist <- dat_sim$bearings[index, ]
    noise_call <- dat_sim$noise_call[index, ]
    noise_random <- dat_sim$noise_random[1:100, ]
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
    f_density <- D ~ 1 #s(dist_to_coast, k = 3, fx = TRUE) # + depth
    
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
      par_det_start <- c(U = log(0.8 / (1 - 0.8)), # should be in (0, 1]
                         B = log(B_loop), # should be in (0, Inf)
                         Q = log(Q_loop - 1)) # should be in (1, Inf)
    } else {
      par_det_start <- c(g0 = log(0.7 / (1 - 0.7)),
                         sigma = log(15000))
    }
    
    # Start values for density function
    par_dens_start <- c("(Intercept)" = -1.5)
    # par_dens_start <- c("(Intercept)" = -2, "dist_to_coast" = 1, "dist_to_coast2" = -1.5)
    # par_dens_start <- c("(Intercept)" = -3, "depth" = 0.5, "depth2" = -2)
    
    # Start values for received level
    par_rl_start <- c(beta_r = log(15),
                      sd_r = log(1))
    
    # Start values for source level
    par_sl_start <- c(mu_s = log(100), # identity
                      sd_s = log(5)) # log
    
    # Start values for bearing
    par_bear_start <- c(kappa = log(10)) # log
    
    par_start <- list(par_det = par_det_start,
                      # par_bear = par_bear_start,
                      par_rl = par_rl_start,
                      par_sl = par_sl_start,
                      par_dens = par_dens_start)
    
    ######################## Run bw_scr() on subset data ##########################
    method = "L-BFGS-B"
    maxit = 100
    TRACE = TRUE
    LSE = TRUE
    
    par <- par_start
    
    #################################bw_ascr step #############################
    ################## Load libraries and perform input checks ###################
    library(dplyr)
    library(mgcv)
    library(matrixStats)
    library(circular)
    
    source("Scripts/Function_v2/HTLikeEstimator.R")
    
    # Input checks
    if (class(dat) != "list") {
      stop()
    }
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
    
    # Check which detection function is to be used, and make sure detection rate
    # or probability at distance = 0 is correctly named. 
    if (dat$det_function == "jano") {
      CORRECT_PAR <- all(c("U", "B", "Q") %in% names(par[["par_det"]]))
      if (!CORRECT_PAR) {
        stop("Incorrect start parameters specified for Janoschek detection function.")
      }
    } else if (dat$det_function == "half-normal") {
      CORRECT_PAR <- all(c("g0", "sigma") %in% names(par[["par_det"]]))
    } else {
      stop(paste0("Detection function specification is ", dat$det_function,
                  ", but should be 'jano' or 'half-normal"))
    }
    
    # Turn par into a named vector with the correct names (optim() requires a vector)
    par <- unlist(par)
    names(par) <- gsub(".*\\.", "", names(par))
    ## TO KEEP SOME PARAMETERS FIXED, USE LINE BELOW TO SELECT WHICH PARAMETERS TO ESTIMATE
    # par <- par[names(par) %in% c("g0", "sigma" , "kappa", "(Intercept)", "dist_to_coast", "dist_to_coast2")] # USE ONLY FOR SIMULATIONS
    # par <- par[names(par) %in% c("U", "B", "Q", "mu_s", "beta_r", "(Intercept)")] # USE ONLY FOR SIMULATIONS
    
    # Create a matrix of distances
    distances <- apply(dat$detectors, 1, function(z) {
      y <- t(dat$cov_density[, c("x", "y")]) - z
      y <- y ^ 2
      return(sqrt(colSums(y)))
    })
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
    
    # not sure if scaling is necessary
    cov_density_scaled <- dat$cov_density
    cov_density_scaled[, -c(1, 2)] <- scale(select(dat$cov_density, -c("x", "y")))
    
    if (USE_BEARINGS) {
      bearings_deg <- circular(dat$bearings, 
                               units = "degrees", 
                               template = "geographics") # 'geographics' means clockwise rotation with 0 at north
      bearings_rad <- conversion.circular(bearings_deg)
      dat[["bearings_rad"]] <- bearings_rad # Add bearings as circular in radians to dat
      
      # Add bearings of all grid points to all detectors to dat
      grid_bearings <- apply(dat$detectors[, c("x", "y")], 1, function(det) {
        x <- dat$cov_density$x - det[1]
        y <- dat$cov_density$y - det[2]
        return(coord2rad(x, y, control.circular = list(template = "geographics")))
      })
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
      
      par <- c(par, gam_par)
    } else {dat[["smooth_terms"]] <- NULL}
    
    dat[["design"]] <- model.matrix(gam_fit)
    dat[["gam_fit"]] <- gam_fit
    
    # Add whether trace is TRUE/FALSE
    dat[["TRACE"]] <- TRACE
    
    if (LSE) {
      fn <- .llkLSEParallelSmooth
    } else {
      stop("Non LSE version not available yet.")
      #fn <- .loglikelihood_with_bearings
    }
    
    # Remove density from par if density is constant STILL NEEDS IMPLEMENTATION
    if (CONSTANT_DENSITY) {
      par <- par[!(names(par) %in% vars)]
    }
    
    
    ######## llkLSEParellelSmooth.R part###########################3333
    llk_evals <- rbind(llk_evals, c(B_loop, Q_loop, .llkLSEParallelSmooth(par = par, dat = dat)))
  }
}
