## =============================================================================
## This script simulates data that can be used to test model accuracy.
## Written to match the style and simplicity of the Rcpp algorithms.
## =============================================================================

## =============================================================================
## 1. LOAD LIBRARIES AND DATA, AND SET CONSTANTS
## -----------------------------------------------------------------------------

## Load the libraries
library(tidyverse)
library(mgcv)
# library(matrixStats)
library(circular)
library(raster)
library(geosphere)
library(Rfast)

## Load the mesh
mesh <- read.csv("Data/alaska_albers_grid_adaptive_levels=2_inner=10k_outer=50k_maxD2C=Inf_area=8450_n=438.csv")
A_x <- mesh$area

## Load detector data
site <- 5
year <- 2010
DASAR <- as.data.frame(read_tsv("Data/DASARs.txt"))
DASAR$long <- -abs(DASAR$long) # Turn longitude from westing to easting 
detectors <- DASAR[DASAR$year == year & DASAR$site == site, c("long", "lat")] 

## Create standardised covariates data.frame
covariates <- dplyr::select(mesh, c(distance_to_coast, depth)) %>% 
  mutate(distance_to_coast2 = distance_to_coast ^ 2, 
         depth = abs(depth),
         depth2 = depth ^ 2,
         logdepth = log(depth)) 
covariates <- scale(covariates, center = FALSE)

## Set constants
FIXED_SL <- FALSE
MIX_BEAR <- TRUE
trunc_level <- 96
n_grid <- nrow(mesh)
n_det <- nrow(detectors)

## =============================================================================

## =============================================================================
## 2. SET PARAMETER VALUES
## -----------------------------------------------------------------------------
## Detection parameter
g0 <- 0.6

## Propagation parameters
beta_r <- 14.5
sd_r <- 4.5

## Source level parameters
mu_s <- 155
if (!FIXED_SL) sd_s <- 5

## Bearing parameters
if (MIX_BEAR) {
  kappa_low <- 0.3
  kappa_high <- 35
  mix_par <- 0.1 # what percentage has low accuracy, ie uses kappa_low?
} else kappa <- 15

## Density parameters
# f_density <- D ~ 1
f_density <- D ~ distance_to_coast + distance_to_coast2
# par_dens <- c("(Intercept" = -1)
par_dens <- c("(Intercept)" = -16,
              "distance_to_coast" = 57,
              "distance_to_coast2" = -68.5)
## =============================================================================

## =============================================================================
## 3. DERIVE INFORMATION FOR SIMULATION
## -----------------------------------------------------------------------------

## Derive distances
distances <- pointDistance(p1 = mesh[, c("long", "lat")], 
                           p2 = detectors[, c("long", "lat")], 
                           lonlat = TRUE)

## Create GAM object
gam_fit <- gam(f_density, data = data.frame(D = 0, covariates))
if (length(gam_fit$smooth) != 0) {
  k <- gam_fit$smooth[[1]]$bs.dim
  gam_par <- rep(0, k - 1)
  smooth_terms <- paste0(gam_fit$smooth[[1]]$label, ".", 1:(k - 1))
  names(gam_par) <- smooth_terms
  
  par <- c(par, gam_par)
} else {smooth_terms <- NULL}

## Derive grid bearings 
# Get degrees bearings from detectors to all grid points 
grid_bearings <- t(apply(mesh[, c("long", "lat")], 1, function(grid_point) {
  bear <- geosphere::bearing(p1 = as.matrix(detectors[, c("long", "lat")]),
                             p2 = grid_point)
}))
# First convert to 'circular' 
grid_bearings <- circular::circular(grid_bearings, template = "geographics", 
                                    modulo = "2pi", units = "degrees")
# Then convert to radians
grid_bearings_rad <- circular::conversion.circular(grid_bearings)

## =============================================================================
## 4. START THE SIMULATION
## -----------------------------------------------------------------------------

## Estimate densities
gam_fit$coefficients <- par_dens
D <- exp(predict(gam_fit, data = model.matrix(gam_fit)))

## Simulate calls per grid cell (assuming Poisson point process)
calls_per_cell <- rpois(n = n_grid, lambda = A_x * D)
n_call <- sum(calls_per_cell)
cat(paste0(n_call, " bowhead whale calls were simulated.\n"))

## Extract the distances and covariates for all of the n simulated calls
call_sample_index <- rep(1:n_grid, times = calls_per_cell)
distances_call <- distances[call_sample_index, ]
all_calls <- mesh[call_sample_index, ]

## Create data objects to store simulated data in
S <- vector(length = n_call)
W <- matrix(nrow = n_call, ncol = n_det)
R <- matrix(nrow = n_call, ncol = n_det)
Y_deg <- matrix(nrow = n_call, ncol = n_det)
Y_rad <- matrix(nrow = n_call, ncol = n_det)

## Simulated characteristics for every call
for (i in 1:n_call) {
  ## Simulate a source level
  S[i] <- ifelse(FIXED_SL, mu_s, rnorm(1, mu_s, sd_s)) 
  
  ## Simulate received levels
  R[i, ] <- S[i] - beta_r * log(distances_call[i, ], base = 10) + 
    rnorm(n_det, 0, sd_r)
  
  ## Simulate detection history
  det_probs <- g0 * (R[i, ] > trunc_level)
  W[i, ] <- det_probs > runif(n_det, 0, 1)
  
  ## Simulate bearings
  bearing <- geosphere::bearing(p1 = as.matrix(detectors[, c("long", "lat")]),
                             p2 = all_calls[i, c("long", "lat")])
  
  bearing_deg <- circular::circular(x = bearing, template = "geographics", 
                                    modulo = "2pi", units = "degrees")
  if (!MIX_BEAR) {
    errors <- circular::rvonmises(n_det,
                                  circular::circular(0, template = "geographics", 
                                                     modulo = "2pi", 
                                                     units = "degrees"), 
                                  kappa, 
                                  list(template = "geographics",  
                                       modulo = "2pi", 
                                       units = "degrees"))
    Y_i <- circular::circular(bearing_deg + errors, template = "geographics", 
                              modulo = "2pi", units = "degrees")
    Y_deg[i, ] <- Y_i
    Y_rad[i, ] <- circular::conversion.circular(Y_i)
  } else {
    ## Get low accuracy errors
    errors_low <- circular::rvonmises(n_det,
                                      circular::circular(0, template = "geographics", 
                                                         modulo = "2pi", 
                                                         units = "degrees"), 
                                      kappa_low, 
                                      list(template = "geographics",  
                                           modulo = "2pi", 
                                           units = "degrees"))
    ## Get high accuracy errors
    errors <- circular::rvonmises(n_det,
                                  circular::circular(0, template = "geographics", 
                                                     modulo = "2pi", 
                                                     units = "degrees"), 
                                  kappa_high, 
                                  list(template = "geographics",  
                                       modulo = "2pi", 
                                       units = "degrees"))
    ## Get indices for low accuracy bearings
    index_low <- runif(n_det, 0, 1) < mix_par
    ## Replace some high accuracy error by low accuracy errors
    errors[index_low] <- errors_low[index_low] 
    
    Y_i <- circular::circular(bearing_deg + errors, template = "geographics", 
                              modulo = "2pi", units = "degrees")
    Y_deg[i, ] <- Y_i
    Y_rad[i, ] <- circular::conversion.circular(Y_i)
  }
  
  ## Only keep information for detectors that were involved in a detection
  R[i, !W[i, ]] <- NA
  Y_deg[i, !W[i, ]] <- NA
  Y_rad[i, !W[i, ]] <- NA
}

## Remove the calls with less than two detections
enough_dets <- Rfast::rowsums(W) > 1
cat(paste(sum(enough_dets), 
          "bowhead whale calls were detected at least twice.\n"))

## =============================================================================
## 5. CREATE LIST WITH ALL SIMULATED DATA AND OTHER INFO
## -----------------------------------------------------------------------------
if (MIX_BEAR) {
  if (FIXED_SL) {
    true_pars = list(mu_s, 
                     g0,
                     sd_r,
                     beta_r,
                     kappa_low,
                     kappa_high,
                     mix_par,
                     par_dens)
  } else {
    true_pars = list(mu_s, 
                     sd_s,
                     g0,
                     sd_r,
                     beta_r,
                     kappa_low,
                     kappa_high,
                     mix_par,
                     par_dens)
  }
} else {
  if (FIXED_SL) {
    true_pars = list(mu_s, 
                     g0,
                     sd_r,
                     beta_r,
                     kappa,
                     par_dens)
  } else {
    true_pars = list(mu_s, 
                     sd_s,
                     g0,
                     sd_r,
                     beta_r,
                     kappa,
                     par_dens)
  }
}

sim_data <- list(det_hist = W[enough_dets, ],
                 bearings_rad = Y_rad[enough_dets, ],
                 received_levels = R[enough_dets, ],
                 grid_bearings_rad = grid_bearings_rad,
                 trunc_level = trunc_level,
                 f_density = f_density,
                 FIXED_SL = FIXED_SL,
                 MIX_BEAR = MIX_BEAR,
                 distances = distances,
                 detectors = detectors,
                 mesh = mesh,
                 A_x = A_x,
                 covariates = covariates,
                 true_pars = true_pars,
                 N_expected = sum(A_x * D),
                 N_realised = n_call)
## =============================================================================

rm(list = setdiff(ls(), c("sim_data", "output", "i")))

