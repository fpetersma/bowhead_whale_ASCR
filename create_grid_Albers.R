# ==============================================================================
# New script for creating a grid with equal spacing and areas based on ideas
# from supervisors on 14-05-2021.
# ==============================================================================

# ==============================================================================
# Load relevant packages
# ==============================================================================

library(oce)        # for coast line info
library(sp)         # for spatial points
library(rgdal)      # for coordinate systems
library(ggplot2)    # for plotting
library(marmap)     # for bathymetry data
library(Rfast)      # for fast matrix operations
library(raster)     # to calculate distances
library(tidyverse)  # for data manipulation

# ==============================================================================
# 
# ==============================================================================

# Alaska Albers Equal Area Conic
# https://spatialreference.org/ref/esri/102006/
#Proj4 string
# +proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154
#   +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs

# Have a quick look
aaeac <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 
+y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
wgs84 <- "+proj=longlat +datum=WGS84 +no_defs"

# data(coastlineWorld)
# mapPlot(coastlineWorld, 
#         longitudelim = c(-150, -136),
#         latitudelim = c(68.5, 72), projection=aaeac, col='gray')

# Get detector locations and convert DASAR
DASAR <- as.data.frame(readr::read_tsv("Data/DASARs.txt")) %>% 
  filter(year == 2014, # correct year
         site == 5, # correct site
         pos != "F") # remove bad detector 
DASAR$long <- -DASAR$long # remove the Westing from longitude units
DASAR_coord <- dplyr::select(DASAR, c(long, lat))  
coordinates(DASAR_coord) <- c("long", "lat")
proj4string(DASAR_coord) <- CRS(wgs84)
alaska_albers <- spTransform(DASAR_coord, CRS=CRS(aaeac))
DASAR_coord <- data.frame("x_aa" = coordinates(alaska_albers)[, 1], 
                          "y_aa" = coordinates(alaska_albers)[, 2],
                          DASAR)

#### OLD METHOD START

# # Base the spacing on distance from centroid. is not perfect, but find for now
# distance <- c(20, 120) # buffer widths from centroid
# area <- c(6.25, 25) # in squared kilometres
# spacing <- sqrt(area * 1000000) # spacing in meters
# centre <- c(x_aa = mean(DASAR_coord$x_aa), # the centroid
#             y_aa = mean(DASAR_coord$y_aa))
# 
# repeats <- distance * 1000 / spacing
# 
# # It is probably easiest to create 2 different grids, a fine and a coarse one, 
# # and then safe all centroids and then stick them together.
# 
# # First create the coarse one using raster
# coarse_lims <- as(extent(centre["x_aa"] - 120 * 1000, 
#                          centre["x_aa"] + 120 * 1000, 
#                          centre["y_aa"] - 120 * 1000, 
#                          centre["y_aa"] + 120 * 1000), "SpatialPolygons")
# ras_coarse <- raster(coarse_lims, resolution = spacing[2])
# # res(ras_coarse)
# # coordinates(ras_coarse)
# plot(coordinates(ras_coarse))
# 
# fine_lims <- as(extent(centre["x_aa"] - distance[1] * 1000, 
#                        centre["x_aa"] + distance[1] * 1000, 
#                        centre["y_aa"] - (distance[1] + 10) * 1000, # add 10 on the y-axis
#                        centre["y_aa"] + (distance[1] + 10) * 1000), # add 10 on the y-axis
#                 "SpatialPolygons")
# ras_fine <- raster(fine_lims, resolution = spacing[1])
# coordinates(ras_fine)
# plot(coordinates(ras_fine))
# 
# # Match the two coordinates
# coords_fine <- data.frame(coordinates(ras_fine), res = spacing[1], area = area[1])
# coords_coarse <- data.frame(coordinates(ras_coarse), res = spacing[2], area = area[2])

### END OLD METHOD START NEW METHOD

distance <- c(20, 100) # buffer widths from centroid
spacing <- c(2.5, 5, 20)
area <- spacing ^ 2 # in squared kilometres (note how they are multiples in spacing)
spacing_m <- spacing * 1000 # spacing in meters
centre <- c(x_aa = mean(DASAR_coord$x_aa), # the centroid
            y_aa = mean(DASAR_coord$y_aa))

# It is probably easiest to create 2 different grids, a fine and a coarse one, 
# and then safe all centroids and then stick them together.

# First create the coarse one using raster
coarse_lims <- as(extent(centre["x_aa"] - distance[2] * 1000, 
                         centre["x_aa"] + distance[2] * 1000, 
                         centre["y_aa"] - distance[2] * 1000, 
                         centre["y_aa"] + distance[2] * 1000), "SpatialPolygons")
ras_coarse <- raster(coarse_lims, resolution = spacing_m[2])
# res(ras_coarse)
# coordinates(ras_coarse)
# plot(coordinates(ras_coarse))

fine_lims <- as(extent(centre["x_aa"] - distance[1] * 1000, 
                       centre["x_aa"] + distance[1] * 1000, 
                       centre["y_aa"] - (distance[1] + 10) * 1000, # add 10 on the y-axis
                       centre["y_aa"] + (distance[1] + 10) * 1000), # add 10 on the y-axis
                "SpatialPolygons")
ras_fine <- raster(fine_lims, resolution = spacing_m[1])
# coordinates(ras_fine)
# plot(coordinates(ras_fine))

# Match the two coordinates
coords_fine <- data.frame(coordinates(ras_fine), res = spacing_m[1], area = area[1])
coords_coarse <- data.frame(coordinates(ras_coarse), res = spacing_m[2], area = area[2])

# ==============================================================================
# Keep a fine grid with a certain buffer of the area, and coarse outside of it 
# ==============================================================================

# Label all coarse grid points with ID and whether they are too far away
coords_coarse$ID <- 1:nrow(coords_coarse)

# Create spatial points object of grid
coarse <- coords_coarse
coordinates(coarse) <- c("x", "y")
proj4string(coarse) <- CRS(aaeac)

# Transform object to WGS84 long lat coordinates
coarse <- spTransform(coarse, CRSobj = CRS(wgs84))
coarse <- as.data.frame(coordinates(coarse)) # extract the coordinates as a data.frame from matrix
colnames(coarse) <- c("long", "lat")

# Distance to the nearest DASAR (m)
distances_DASAR <- pointDistance(p1 = coarse[, c("long", "lat")], 
                                 p2 = DASAR[, c("long", "lat")], 
                                 lonlat = TRUE)
coords_coarse$to_nearest_DASAR <- rowMins(distances_DASAR, value = TRUE)
coords_coarse$inner_mesh <- coords_coarse$to_nearest_DASAR < 10000

rm(coarse, distances_DASAR)

# Create column for coarse ID in coords_fine
coords_fine$coarse_ID <- NA

# For every fine grid point, find which coarse grid point it falls in
for (j in 1:nrow(coords_coarse)) {
  xy_coarse <- coords_coarse[j, ]
  x_bounds <- c(xy_coarse$x - xy_coarse$res, xy_coarse$x + xy_coarse$res)
  y_bounds <- c(xy_coarse$y - xy_coarse$res, xy_coarse$y + xy_coarse$res)
  
  for (i in 1:nrow(coords_fine)) {
    if (coords_fine[i, "x"] >= x_bounds[1] & 
        coords_fine[i, "x"] <= x_bounds[2] & 
        coords_fine[i, "y"] >= y_bounds[1] & 
        coords_fine[i, "y"] <= y_bounds[2]) {
      coords_fine[i, c("coarse_ID", "inner_mesh")] <- c(xy_coarse$ID, 
                                                        xy_coarse$inner_mesh) 
      
    }
  }
}

# ggplot(coords_coarse) + geom_point(aes(x = x, y = y, colour = inner_mesh))
# ggplot(coords_fine) + geom_point(aes(x = x, y =y, colour = inner_mesh))
# 
# plot(coords_coarse[, c("x","y")])
# points(coords_fine[coords_fine$inner_mesh == TRUE, c("x","y")])

coords_fine$x <- coords_fine$x + spacing[1] * 1000
coords_fine$y <- coords_fine$y - spacing[1] * 1000

grid <- rbind(coords_coarse[coords_coarse$inner_mesh == FALSE, c("x", "y", "res", "area")], 
              coords_fine[coords_fine$inner_mesh == TRUE, c("x", "y", "res", "area")])


plot(grid[, 1:2], col = grid$area)
points(DASAR_coord[, 1:2], col = "red", pch = 20)

### OLD STUFF
# xlims <- c(seq(from = centre$x_aa - distance[2], 
#                to = centre$x_aa - distance[1], 
#                by = spacing[2]),
#            seq(from = centre$x_aa - distance[1], 
#                to = centre$x_aa + distance[1], 
#                by = spacing[1]),
#            seq(from = centre$x_aa + distance[1], 
#                to = centre$x_aa + distance[2], 
#                by = spacing[2]))
# 
# xlims <- seq(from = centre$x_aa - distance * 1000,
#              to = centre$x_aa + distance * 1000,
#              by = spacing)
# ylims <- seq(from = centre$y_aa - distance * 1000,
#              to = centre$y_aa + distance * 1000,
#              by = spacing)
# grid <- expand.grid(x_aa = xlims, y_aa = ylims)
# ggplot(data = grid) + 
#   geom_point(mapping = aes(x = x_aa, y = y_aa)) + 
#   coord_equal()

# ==============================================================================
# Once grid has been created, convert to longlat coords
# ==============================================================================

# Create spatial points object of grid
grid_aa <- grid
coordinates(grid_aa) <- c("x", "y")
proj4string(grid_aa) <- CRS(aaeac)

# Transform object to WGS84 long lat coordinates
grid_ll <- spTransform(grid_aa, CRSobj = CRS(wgs84))
mesh <- as.data.frame(coordinates(grid_ll)) # extract the coordinates as a data.frame from matrix
colnames(mesh) <- c("long", "lat")
# mesh$long <- -mesh$long # remove the West (W) from longitude

# ==============================================================================
# Get bathymetry info for every gridpoint
# ==============================================================================

# Get bathymetry data at highest resolution
long_lim <- c(-150, -136) # The entire width of UTM zone 6 and 7, roughly 6 degrees
# east and west of the array, plus a bit more in the east
lat_lim <- c(68.5, 72) # Should be enough latitude

# Extract bathymetry data from NOAA
bathy <- getNOAA.bathy(lon1 = long_lim[1], lon2 = long_lim[2], 
                       lat1 = lat_lim[1], lat2 = lat_lim[2],
                       resolution = 1, keep = TRUE) # res = 1/30
mesh <- get.depth(bathy, x = mesh$long, y = mesh$lat, locator = FALSE)
colnames(mesh)[1] <- "long"

# ==============================================================================
# Get all distances and add back area and spacing
# ==============================================================================

mesh <- cbind(mesh, res = grid$res, area = grid$area)

# To the nearest DASAR (m)
distances_DASAR <- pointDistance(p1 = mesh[, c("long", "lat")], 
                                 p2 = DASAR[, c("long", "lat")], 
                                 lonlat = TRUE)
mesh$to_nearest_DASAR <- rowMins(distances_DASAR, value = TRUE)

# To the coast (m)
# Get coast line data
data(coastlineWorldFine, package = "ocedata")
coast <- na.omit(as.data.frame(coastlineWorldFine@data)) %>% 
  filter(longitude >= long_lim[1],
         longitude <= long_lim[2],
         latitude >= lat_lim[1],
         latitude <= lat_lim[2])
colnames(coast) <- c("long", "lat")

# Derive shortest distance to coast
all_distances_to_coast <- pointDistance(p1 = mesh[, c("long", "lat")], 
                                        p2 = coast, 
                                        lonlat = TRUE)
mesh$distance_to_coast <- rowMins(all_distances_to_coast, value = TRUE)
# ==============================================================================
# Remove parts of mesh that shouldn't be evaluated
# ==============================================================================

# add back Alaska Albers projection coordinates
mesh$x_aa <- grid$x
mesh$y_aa <- grid$y

mesh_filtered <- mesh %>% 
  filter(depth < 0, # Remove coordinates on land
         to_nearest_DASAR <= 50000) # Remove coordinates that are too far away

# ==============================================================================
# Evaluate current mesh
# ==============================================================================
# Visually check the new adaptive mesh
p_depth <- ggplot(data = mesh_filtered) +
  geom_point(mapping = aes(x = x_aa, y = y_aa, colour = depth)) + 
  labs(colour = "Depth") + coord_equal() + 
  geom_point(data = DASAR_coord, mapping = aes(x = x_aa, y = y_aa), colour = "red")
# + geom_line(data = coast_x_aay_aa, mapping = aes(x = x_aa, y = y_aa, size = "Coastline (ocedata)"), 
#           linetype = 1) + scale_size_manual(values = 1)
p_distance_coast<- ggplot(data = mesh_filtered) +
  geom_point(mapping = aes(x = x_aa, y = y_aa, colour = distance_to_coast)) + 
  labs(colour = "D2Coast") + coord_equal() + 
  geom_point(data = DASAR_coord, mapping = aes(x = x_aa, y = y_aa), colour = "red")
# + geom_line(data = coast_x_aay_aa, mapping = aes(x = x_aa, y = y_aa, size = "Coastline (ocedata)"), 
#           linetype = 1) + scale_size_manual(values = 1)
p_distance_DASAR<- ggplot(data = mesh_filtered) +
  geom_point(mapping = aes(x = x_aa, y = y_aa, colour = to_nearest_DASAR)) + 
  labs(colour = "D2DASAR") + coord_equal() + 
  geom_point(data = DASAR_coord, mapping = aes(x = x_aa, y = y_aa), colour = "red")
# + geom_line(data = coast_x_aay_aa, mapping = aes(x = x_aa, y = y_aa, size = "Coastline (ocedata)"), 
#           linetype = 1) + scale_size_manual(values = 1)
p_area <- ggplot(data = mesh_filtered) +
  geom_point(mapping = aes(x = x_aa, y = y_aa, colour = area)) + 
  labs(colour = "Area") + coord_equal() + 
  geom_point(data = DASAR_coord, mapping = aes(x = x_aa, y = y_aa), colour = "red")
# + geom_line(data = coast_x_aay_aa, mapping = aes(x = x_aa, y = y_aa, size = "Coastline (ocedata)"), 
#           linetype = 1) + scale_size_manual(values = 1)
gridExtra::grid.arrange(p_depth, p_distance_coast, p_distance_DASAR, p_area,
                        ncol = 1, nrow = 4, 
                        top = "The adaptive grid with four different colourings \n(depth, distance_to_coast, distance_to_DASAR, and area)\n ")

# mapPlot(mesh_filtered$long, latitude = mesh_filtered$lat, projection = aaeac, type = "p")
mapPlot(coastlineWorldFine,
        longitudelim = c(-147, -140),
        latitudelim = c(69.5, 71.5), projection=aaeac, col='gray')

points(mesh_filtered$x_aa, mesh_filtered$y_aa)
points(DASAR_coord$x_aa, DASAR_coord$y_aa, col = "red", pch = 20)

# Write .csv file =========================================================  ##
write.csv(x = mesh_filtered[, c("long", "lat", "area", "depth", "distance_to_coast", "to_nearest_DASAR")],
          file = "Data/alaska_albers_grid_adaptive_levels=2_inner=10k_outer=50k_maxD2C=Inf_area=8450_n=438.csv",
          row.names = FALSE)
