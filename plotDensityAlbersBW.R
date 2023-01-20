## IMPORTANT https://r-spatial.org/r/2018/10/25/ggplot2-sf.html

## This is a black and white version of plotDensityAlbers.R, and probably better

## Load files to try things
mesh <- read.csv("Data/alaska_albers_grid_adaptive_levels=2_inner=10k_outer=50k_maxD2C=Inf_area=8450_n=438.csv")
load(file =  "Real data output/fits_1_35_nlminb_n=443.RData" )
# d <- fit_obj_varDens$D

plotDensityAlbers <- function(d, mesh, name = "density_plot") {
  # ==============================================================================
  # Load relevant packages
  # ==============================================================================
  
  library(oce)        # for coast line info
  library(sp)         # for spatial points
  library(rgdal)      # for coordinate systems
  library(ggplot2)    # for plotting
  # library(marmap)     # for bathymetry data
  # library(ggpattern)  # for patterns as fillers
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
                                                           colour = density), 
               size = 5.5,
               # alpha = 0.7,
               shape = 15) + 
    geom_point(data = mesh[d$area == 6.25, ],  mapping = aes(x = long.1, 
                                                             y = lat.1, 
                                                             colour = density), 
               size = 2.6,
               # alpha = 0.7,
               shape = 15) + 
    # scale_colour_distiller(type = "seq",
    #                       direction = 1,
    #                       palette = "Greys") +
    # guides(alpha = "none") +
    # scale_colour_gradientn(colours =c("grey85", "grey85", "grey30", "grey0"),values = c(0, 0.01, 0.1, 1)) + 
    scale_colour_gradient(low = alpha("grey90", 0.8),  
                          high = alpha("grey0", 0.8), 
                          trans= "sqrt", 
                          name = "Calls/km2") +
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
    theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", 
                                          linewidth = 0.5), 
          panel.background = element_rect(fill = "white")) +
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
    # labs(colour = "Calls/km2") + 
    geom_point(data = DASAR, 
               mapping = aes(x = long.1, y = lat.1),
               shape = 17, 
               size = 2, show.legend = TRUE) + 
    # ggtitle("Map of study area, with the density surface overlayed") +
    xlab("Longitude") + 
    ylab("Latitude") #+ geom_point(data=mesh[c(290, 275), ], 
                                  # mapping=aes(x = long.1, y = lat.1),
                                  # colour = 'red') + geom_point(data=mesh[c(277,278,279,289, 291), ], 
                                  #                              mapping=aes(x = long.1, y = lat.1),
                                  #                              colour = 'blue')
                                  # 
  
  return(p)
}

for (ii in 33) {
  p <- plotDensityAlbers(d = fits[[ii]]$D, mesh = mesh)
  
  ggsave(filename = paste0("grey_map_", ii, ".pdf"), plot = p,  height = 12, width = 8.39)
  ggsave(filename = paste0("grey_map_", ii, ".png"), plot = p, width = 8.275, height = 7, dpi = "retina") # these dimensions work!
  
}
p

# save(list = "p", file = "density_plot.RData")


plotDensityAlbers(d = d, mesh = mesh)

ggsave("map_31.pdf", height = 12, width = 8.39)
ggsave("map_31.png", width = 8.275, height = 12, dpi = "retina") # these dimensions work!
