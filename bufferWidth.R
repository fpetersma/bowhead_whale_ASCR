bufferWidth <- function(p,
                        trunc_level,
                        g0, 
                        beta_r, 
                        sd_r, 
                        mu_s, 
                        sd_s = NULL, 
                        noise = NULL,
                        SINGLE_SL = TRUE, 
                        WITH_NOISE = FALSE) {
  # Get the correct quantile for the normal cdf, given p
  quant <- qnorm(1 - p)
  
  # Derive the distance at which that p is achieved
  distance <- 10 ^ ((quant * sd_r - trunc_level + mu_s) / beta_r)
  
  return(distance)
  # Assuming all detectors are at this distance, what is p.? 
  # RECALL: we are only interested in calls that were detected at least twice!)
  
}

bufferWidth(p = 0.01,
            trunc_level = 96,
            g0 = 0.62,
            beta_r = 17.36,
            sd_r = 4.85,
            mu_s = 165)

# ==============================================================================
# Create a function that calculates p. for every grid point of the mesh
# ------------------------------------------------------------------------------
pDetected <- function(mesh,
                      detectors,
                      trunc_level,
                      min_no_detections = 2, 
                      g0, 
                      beta_r, 
                      sd_r, 
                      mu_s, 
                      sd_s = NULL, 
                      noise = NULL,
                      SINGLE_SL = TRUE, 
                      WITH_NOISE = FALSE) {
  
  source("Scripts/bowhead whales/hidden_functions.R")
  
  # Create a matrix of distances from mesh grid points to detectors
  distances <- pointDistance(p1 = mesh[, c("long", "lat")], 
                             p2 = detectors[, c("long", "lat")], 
                             lonlat = TRUE)
  
  
  
  # calculate the detection probabilities
  if (SINGLE_SL & !WITH_NOISE) {
    # calculate the received levels
    levels <- mu_s - beta_r * log(distances, base = 10)
    
    det_probs <- g0 * (1 - pnorm((trunc_level - levels) / sd_r))
    
    p. <- .detected(probs = det_probs, 
                    min_no_detections = min_no_detections)
  }
  if (!SINGLE_SL & !WITH_NOISE) {
    p. <- rowSums(sapply(100:200, function(x) {
      levels <- x - beta_r * log(distances, base = 10)
      det_probs <- g0 * (1 - pnorm((trunc_level - levels) / sd_r))
      p._temp <- .detected(probs = det_probs, 
                      min_no_detections = min_no_detections)
      p <- dnorm(x, mu_s, sd_s)
      
      return(p._temp * p)
    }))
  }
  
  return(p.)
  
}

# load libraries
library(readr)
library(dplyr)
library(mgcv)
library(matrixStats)
library(circular)
library(raster)

# Get the mesh
mesh_file <- "Data/grid_adaptive_levels=2_bounds=10k_maxD2C=Inf_maxD2A=60k_area=11461.8_n=737+567=1304.csv"
mesh <- read.csv(mesh_file)

# Get the detectors
DASAR <- as.data.frame(read_tsv("Data/DASARs.txt"))
DASAR$long <- -abs(DASAR$long) # Turn longitude from westing to easting (more common)
detectors <- DASAR[DASAR$year == 2010 & DASAR$site == 5, c("long", "lat")]

output <- pDetected(mesh = mesh, 
                    detectors = detectors, 
                    trunc_level = 96,
                    min_no_detections = 2,
                    g0 = 0.505,
                    beta_r = 17.36,
                    sd_r = 2.36,
                    mu_s = 158.5,
                    sd_s = 5.23,
                    SINGLE_SL = FALSE)

probs <- data.frame(det_prob = output, mesh)


library(ggplot2)
fig <- ggplot(data = probs) +
  geom_point(mapping = aes(x = long, y = lat, colour = det_prob >= 0.001), 
             alpha = 0.5, shape = 16) + 
  labs(colour = "Within 0.1% detection") +
  coord_equal() + 
  # scale_colour_manual(values = c("darkred", "darkgreen")) +
  geom_point(data = detectors, 
             mapping = aes(x = long, y = lat, colour = "red"))
fig

points(x = detectors$long, y = detectors$lat, col = "red")

#### DOES NOT WORK FROM HERE


library(plotly)
library(reshape2)


casted <- dcast(probs, long ~ lat, value.var = "det_prob") # does not work

fig <- plot_ly(z = probs$det_prob, 
               y = probs$long, 
               x = probs$lat)
fig <- fig %>% add_surface()
fig

library(plot3D)
surf3D(z = probs$det_prob, 
       y = probs$long, 
       x = probs$lat)
