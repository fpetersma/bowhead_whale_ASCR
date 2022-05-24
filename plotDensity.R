## IMPORTANT https://r-spatial.org/r/2018/10/25/ggplot2-sf.html

## Load files to try things
mesh <- read.csv("Data/alaska_albers_grid_adaptive_levels=2_inner=10k_outer=50k_maxD2C=Inf_area=8450_n=438.csv")
# load(file =  "Simulation output/Exploratory fits to set simulation parameters/fit_obj_varDens_center=F.RData" )
# d <- fit_obj_varDens$D
d <- fits[[31]]$D

plotDensityAlbers <- function(d, mesh, name = "density_plot") {
  # ==============================================================================
  # Load relevant packages
  # ==============================================================================
  
  library(oce)        # for coast line info
  library(sp)         # for spatial points
  library(rgdal)      # for coordinate systems
  library(ggplot2)    # for plotting
  # library(marmap)     # for bathymetry data
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
  
  # Add density to mesh file
  mesh$density <- d$density

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
  
  range(mesh$long.1)
  range(mesh$lat.1)
  ## ===========================================================================
  ## Print current mesh
  ## ===========================================================================
  
  
  theme_set(theme_bw())
  
  
  world <- ne_countries(scale = "large", returnclass = "sf")
  class(world)
  
  p <- ggplot(data = world) +
    # First plot the points
    geom_point(data = mesh[d$area == 25, ],  mapping = aes(x = long.1, 
                                                           y = lat.1, 
                                                           colour = density,
                                                           alpha = -to_nearest_DASAR), 
               size = 5.5,
               shape = 15) + 
    geom_point(data = mesh[d$area == 6.25, ],  mapping = aes(x = long.1, 
                                                             y = lat.1, 
                                                             colour = density), 
               size = 2.6,
               shape = 15) + 
    scale_colour_viridis_c() +
    guides(alpha = FALSE) +
    # scale_colour_gradient(low = "white", high = "darkgrey") + 
    # scale_fill_viridis_c(option = "plasma", trans = "sqrt") +
  
    # The look of the map
    geom_sf(fill= "antiquewhite", colour = "black") + 
    annotation_scale(location = "bl", width_hint = 0.5) +
    annotation_north_arrow(location = "bl", 
                           which_north = "true", 
                           pad_x = unit(0.75, "in"), 
                           pad_y = unit(0.5, "in"), 
                           style = north_arrow_fancy_orienteering) + 
    coord_sf(xlim = c(320000, 490000), ylim = c(2220000, 2380000), crs = aaeac) +
    theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), 
          panel.background = element_rect(fill = "aliceblue")) +
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
    labs(colour = "Density\n(calls/km2)") + 
    geom_point(data = DASAR, 
               mapping = aes(x = long.1, y = lat.1), 
               colour = "red",
               shape = 17, 
               size = 3, show.legend = TRUE) + 
    xlab("Longitude") + 
    ylab("Latitude") +
    ggtitle("Map of study area, with the density surface overlayed") 
  
  return(p)
  
  # ## Visually check the new adaptive mesh
  # p_density <- ggplot() +
  #   geom_point(data = mesh[d$area == 25, ],  mapping = aes(x = long.1, y = lat.1, colour = density), 
  #              size = 6.5,
  #              shape = 15) + 
  #   geom_point(data = mesh[d$area == 6.25, ],  mapping = aes(x = long.1, y = lat.1, colour = density), 
  #              size = 3,
  #              shape = 15) + 
  #   labs(colour = "Density") + coord_equal() + 
  #   geom_point(data = DASAR, 
  #              mapping = aes(x = long.1, y = lat.1), 
  #              colour = "red",
  #              shape = 17, size = 3)
  # p_density
  # 
  # # To the coast (m)
  # # Get coast line data
  # data(coastlineWorldFine, package = "ocedata")
  # 
  # mapPlot(coastlineWorldFine,
  #         longitudelim = c(-147, -139),
  #         latitudelim = c(69.5, 71.5), projection=aaeac, col='gray')
  # 
  # points(mesh$long.1, mesh$lat.1, col = "red")
  # points(DASAR_coord$x_aa, DASAR_coord$y_aa, col = "red", pch = 20)
  # 
  # 
  # 
  # p <- ggplot()
  # # p <- p + geom_raster(data = bathy, aes(x = x, y = y, fill = depth_bins), interpolate = TRUE, alpha = 0.75)
  # p <- p + geom_map(data = as.data.frame(coastlineWorldFine@data), aes(x = longitude, y = latitude))
  # p <- p + ylim(c(69.5, 71.5)) + xlim(c(-147, -139))
  # p <- p + coord_proj(wgs84)
  # p
}

for (ii in 1:29) {
  p <- plotDensityAlbers(d = fits[[ii]]$D, mesh = mesh)
  
  ggsave(filename = paste0("map_", ii, ".pdf"), plot = p,  height = 12, width = 8.39)
  ggsave(filename = paste0("map_", ii, ".png"), plot = p, width = 8.275, height = 12, dpi = "retina") # these dimensions work!
  
}


plotDensityAlbers(d = d, mesh = mesh)

ggsave("map_31.pdf", height = 12, width = 8.39)
ggsave("map_31.png", width = 8.275, height = 12, dpi = "retina") # these dimensions work!

# 
# plotDensity <- function(d) {
#   # Takes either a matrix of densities with x and y coordinates, or an
#   # object of class "bw_ascr_fit" containing a density matrix, and creates a 
#   # density/heat map. 
#   
#   # Load necessary library 'lattice'
#   library(lattice)
#   library(raster)
#   library(ggplot2)
#   
#   if (class(d) == "bwASCR_model") {
#     d <- d$density
#   }
# 
#   # ggplot(data = d, mapping = aes(x = long, y = lat, fill = density)) +
#   #   geom_tile()
#   
#   # Create the heat map
#   par(mar=c(3,4,2,2))
#   levelplot(density ~ long * lat,
#             data = d,
#             pretty = TRUE,
#             xlab = "Longitude easting",
#             ylab = "Latitude northing",
#             col.regions = heat.colors(100)[length(heat.colors(100)):1],
#             main = "")
# 
# }
