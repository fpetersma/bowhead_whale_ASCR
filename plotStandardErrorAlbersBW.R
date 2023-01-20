## ========================================================================== ##
## create_SE_density_plot.R                                                   ##
##                                                                            ##
## This script creates the standard error map from the best model             ##
## ========================================================================== ##

## Load files to try things
mesh <- read.csv("Data/alaska_albers_grid_adaptive_levels=2_inner=10k_outer=50k_maxD2C=Inf_area=8450_n=438.csv")

## Load bootstraps
load("Real data output/results_model_33_1_999.RData")

## Load best model
# load("Real data output/fits_1_21_nlminb_reltol=1e-8_n=443.RData")
load("Real data output/fits_1_35_nlminb_n=443.RData")

## Combine and clean
results <- c(results, fits[33])
rm(fits)

## EVERY FIT CONTAINS A DENSITY MAP, SO EXTRACT THE STANDARD ERROR FOR EVERY ONE
densities <- sapply(results, function(x) x$D$density)

## Standard error
## Absolute measure of spread, sensitive to outliers
se_densities <- data.frame(area = results[[1]]$D$area,
                           spread = apply(densities, 1, sd))

## Coefficient of variation
## Is a relative measure of spread, sensitive to outliers
cv_densities <- data.frame(area = results[[1]]$D$area,
                           spread = apply(densities, 1, function(r) {
                             sd(r) / mean(r) * 100
                           }))

iqr_densities <- data.frame(area = results[[1]]$D$area,
                           spread = apply(densities, 1, function(r) {
                             qt <- quantile(r)
                             iqr <- (qt[4] - qt[2])
                           }))

## Quartile coefficient of dispersion (Q3 - Q1) / (Q1 + Q3)
## Is a relative measure of spread unaffected by outliers
qcd_densities <- data.frame(area = results[[1]]$D$area,
                            spread = apply(densities, 1, function(r) {
                              qt <- quantile(r)
                              qcd <- (qt[4] - qt[2]) / (qt[4] + qt[2])
                            }))
# turn NA values to 0
qcd_densities$spread[is.nan(qcd_densities$spread)] <- 0

## Create function similar to the one in plotDensityAlbersBW.R
plotSE <- function(spread, mesh, legend_title) {
  # ==============================================================================
  # Load relevant packages
  # ==============================================================================
  
  library(oce)        # for coast line info
  library(sp)         # for spatial points
  library(rgdal)      # for coordinate systems
  library(ggplot2)    # for plotting
  library(Rfast)      # for fast matrix operations
  library(raster)     # to calculate distances
  library(tidyverse)  # for data manipulation
  library(ggalt)
  
  library("sf")
  library("rgeos")
  library("ggspatial")
  
  library("rnaturalearth")
  library("rnaturalearthdata")
  library("rnaturalearthhires") 
  
  ## Create relevant projections
  aaeac <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 
+y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
  wgs84 <- "+proj=longlat +datum=WGS84 +no_defs"
  
  
  # Get detector locations and convert DASAR
  DASAR <- as.data.frame(readr::read_tsv("Data/DASARs.txt")) %>% 
    filter(year == 2014, # correct year
           site == 5, # correct site
           pos != "F") # remove bad detector 
  DASAR$long <- -DASAR$long # remove the Westing from longitude units
  DASAR_coord <- dplyr::select(DASAR, c(long, lat))  
  coordinates(DASAR_coord) <- c("long", "lat")
  proj4string(DASAR_coord) <- CRS(wgs84)
  DASAR_coord_albers <- spTransform(DASAR_coord, CRS=CRS(aaeac))
  DASAR <- data.frame(DASAR, DASAR_coord_albers)
  
  # Add spread to mesh file
  mesh$spread <- spread$spread
  
  ## Convert coordinates to Albers 
  mesh_coord <- dplyr::select(mesh, c(long, lat))  
  coordinates(mesh_coord) <- c("long", "lat")
  proj4string(mesh_coord) <- CRS(wgs84)
  mesh_coord_albers <- spTransform(mesh_coord, CRS=CRS(aaeac))
  
  ## Add coordinates to original mesh
  mesh <- data.frame(mesh, mesh_coord_albers)
  
  # Distance to the nearest DASAR (m)
  distances_DASAR <- pointDistance(p1 = mesh[, c("long", "lat")], 
                                   p2 = DASAR[, c("long", "lat")], 
                                   lonlat = TRUE)
  mesh$to_nearest_DASAR <- rowMins(distances_DASAR, value = TRUE)
  
  # range(mesh$long.1)
  # range(mesh$lat.1)
  ## ===========================================================================
  ## Print current mesh
  ## ===========================================================================
  theme_set(theme_bw())
  
  world <- ne_countries(scale = "large", returnclass = "sf")
  class(world)
  
  p <- ggplot(data = world) +
    # First plot the points
    geom_point(data = mesh[spread$area == 25, ],  mapping = aes(x = long.1,
                                                           y = lat.1,
                                                           colour = spread),
               size = 5.5,
               # alpha = 0.7,
               shape = 15) +
    geom_point(data = mesh[spread$area == 6.25, ],  mapping = aes(x = long.1,
                                                             y = lat.1,
                                                             colour = spread),
               size = 2.6,
               # alpha = 0.7,
               shape = 15) +
    # guides(alpha = "none") +
    # scale_colour_gradient(low = alpha("grey85", 0.8), 
    #                       high = alpha("grey0", 0.8), # Define alpha here instead of above to have legend match it
    #                       name = legend_title) +
    scale_colour_gradient(low = alpha("grey90", 0.8),  
                          high = alpha("grey0", 0.8),  
                          name = legend_title) +
    # scale_fill_viridis_c(option = "plasma", trans = "sqrt") +
    
    # The look of the map
    geom_sf(fill= "white", colour = "black") + 
    annotation_scale(location = "bl", width_hint = 0.5) +
    annotation_north_arrow(location = "bl", 
                           which_north = "true", 
                           pad_x = unit(0.75, "in"), 
                           pad_y = unit(0.5, "in"), 
                           style = north_arrow_fancy_orienteering) + 
    coord_sf(xlim = c(320000, 490000), ylim = c(2220000, 2380000), crs = aaeac) +
    theme(panel.grid.major = element_line(color = gray(.5), 
                                          linetype = "dashed", 
                                          linewidth = 0.5), 
          panel.background = element_rect(fill = "white")) +
          # panel.background = element_rect(fill = NA),
          # panel.ontop = TRUE) +
    annotate(geom = "text", x = 350000, y = 2350000, 
             label = "Beaufort Sea", 
             fontface = "italic", 
             color = "grey22", 
             size = 6) +
    annotate(geom = "text", x = 420000, y = 2240000, 
             label = "Canada", 
             fontface = "italic", 
             color = "grey22", 
             size = 6) +
    # labs(colour = legend_title) + 
    geom_point(data = DASAR, 
               mapping = aes(x = long.1, y = lat.1),
               shape = 17, 
               size = 2, show.legend = TRUE) + 
    # ggtitle("Map of study area, with the density surface overlayed") +
    xlab("Longitude") + 
    ylab("Latitude")
  # p$layers <- c(geom_point(data = mesh[spread$area == 25, ],  mapping = aes(x = long.1,
  #                                                                                         y = lat.1,
  #                                                                                         colour = density),
  #                                        size = 5.5,
  #                                        alpha = 0.7,
  #                                        shape = 15) +
  #                            geom_point(data = mesh[spread$area == 6.25, ],  mapping = aes(x = long.1,
  #                                                                                          y = lat.1,
  #                                                                                          colour = density),
  #                                       size = 2.6,
  #                                       alpha = 0.7,
  #                                       shape = 15),
  #                          p$layers)

  
  return(p)
}

## Print the plot
plotSE(spread = se_densities, mesh = mesh, legend_title = "SE")
ggsave("SE_map_8.eps", width = 7.9, height = 12, dpi = "retina") # these dimensions work!

plotSE(spread = cv_densities, mesh = mesh, legend_title = "CV (%)")
ggsave("CV_map_8.png", width = 8.0, height = 12, dpi = "retina") # these dimensions work!

plotSE(spread = iqr_densities, mesh = mesh, legend_title = "IQR")
ggsave("IQR_map_8.png", width = 8.1, height = 12, dpi = "retina") # these dimensions work!

qcd_plot <- plotSE(spread = qcd_densities, mesh = mesh, legend_title = "QCD")
ggsave("QCD_map_8.png", width = 8.1, height = 12, dpi = "retina") # these dimensions work!

## Make one figure
load("C:/Users/felix/OneDrive - University of St Andrews/Documents/University of St Andrews/PhD/Bowhead Whales/JABES paper/qcd_plot_8.RData")
load("C:/Users/felix/OneDrive - University of St Andrews/Documents/University of St Andrews/PhD/Bowhead Whales/JABES paper/density_plot_8.RData")

two_maps <- gridExtra::grid.arrange(p, qcd_plot, ncol = 2)
ggsave(plot = two_maps, filename = "images/two_maps_qcd_33.png", width = 16.15, height = 10, dpi = "retina") # these dimensions work!
