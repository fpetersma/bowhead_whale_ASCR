##  ========================================================================  ##
##  Create a grid with two resolutions: fine and coarse.                      ##
##                                                                            ##
##  Resolutions: 1 = 1/60, 2 = 1/30, 3 = 1/20, 4 = 7/104 = 1/15 (roughly)     ##
##               5 = 5/60, 6 = 6/60, 7 = 7/60, etc etc
##  ========================================================================  ##

##  Load libraries =========================================================  ##
library(raster)
library(geosphere)
library(maptools)
library(marmap)
library(readr)
library(oce)
library(rgdal)
library(circular)
library(tidyverse)
library(reshape2)

DASAR_bound <- 100000
coast_bound <- 200000

##  Create the map =============================================================

# Load the DASAR locations
DASAR <- as.data.frame(read_tsv("Data/DASARs.txt"))
detectors <-  DASAR %>% 
  filter(year == 2010, site == 5) %>% 
  select(long, lat) %>%
  transmute(long = -long, lat = lat)

x_rangeDASAR <- range(detectors$long)
y_rangeDASAR <- range(detectors$lat)

# Site 5 is in zone 7, site 3 is in zone 6
lon <- c(-150, -136) # The entire width of UTM zone 6 and 7, roughly 6 degrees
# east and west of the array, plus a bit more in the east
lat <- c(68.5, 72) # Hopefully enough

# Extract bathymetry data from NOAA
bathy <- getNOAA.bathy(lon1 = lon[1], lon2 = lon[2], 
                       lat1 = lat[1], lat2 = lat[2],
                       resolution = 3, keep = TRUE)
raster_longlat <- as.raster(bathy)

# Proceed to converting to matrices
coord <- coordinates(raster_longlat)
matrix_longlat <- as.matrix(raster_longlat)
rownames(matrix_longlat) <- rev(sort(unique(coord[, 2])))
colnames(matrix_longlat) <- sort(unique(coord[, 1]))

# Set all land to NA
matrix_longlat[matrix_longlat >= 0] <- NA

## Get the coastlines using oce ============================================  ##
data(coastlineWorldFine, package = "ocedata")
coast <- na.omit(as.data.frame(coastlineWorldFine@data))
coast_longlat <- coast[coast$longitude >= lon[1] &
                         coast$longitude <= lon[2] & 
                         coast$latitude >= lat[1] &
                         coast$latitude <= lat[2], ]
colnames(coast_longlat) <- c("long", "lat")

#--------------------# Calculate the distances #-------------------------------#
melt_longlat <- melt(matrix_longlat, varnames = c("lat", "long"), 
                     value.name = "depth")
melt_longlat <- na.omit(melt_longlat)[, c("long", "lat", "depth")]
distances_longlat <- pointDistance(p1 = melt_longlat[, c("long", "lat")], 
                                   p2 = detectors[, c("long", "lat")], 
                                   lonlat = TRUE)
melt_longlat$to_nearest_DASAR <- Rfast::rowMins(distances_longlat, value = TRUE)

#-----------------# Calculate the distance to coast #--------------------------#
# My own code, as I am not sure if I trust dist2Line()
all_distances_to_coast <- pointDistance(p1 = melt_longlat[, c("long", "lat")], 
                                        p2 = coast_longlat, 
                                        lonlat = TRUE)
distance_to_coast <- Rfast::rowMins(all_distances_to_coast, value = TRUE)

melt_longlat <- cbind(melt_longlat, distance_to_coast = distance_to_coast)

mesh <- subset(melt_longlat, subset = to_nearest_DASAR <= DASAR_bound & 
                 distance_to_coast <= coast_bound)

# Transform the circles to rasters to get the area for every cell. As the earth
# is not flat (duh), the area per cell varies not just between circles, but also
# between different latitude (and longitudes maybe, but I don't think so). Either 
# way, we need an individual area measurement for every grid cells we integrate 
# over. 
raster_1 <- rasterFromXYZ(mesh[, 1:3], 
                          crs = "+proj=longlat +datum=WGS84 +no_defs")

# Here we recreate the melted data we had before. 
df_1 <- na.omit(cbind(coordinates(raster_1), 
                      depth = values(raster_1), 
                      area = area(raster_1)[]))


## Add the two meshes together ============================================== ##
mesh_adaptive <- as.data.frame(df_1)
colnames(mesh_adaptive) <- c("long", "lat", "depth", "area")
sum(mesh_adaptive$area) 

# Add distance to coast and distance to DASAR, as these are now gone.
# Distance to DASAR
distances_adaptive <- pointDistance(p1 = mesh_adaptive[, c("long", "lat")], 
                                    p2 = detectors[, c("long", "lat")], 
                                    lonlat = TRUE)
mesh_adaptive$to_nearest_DASAR <- Rfast::rowMins(distances_adaptive, value = TRUE)

# Distance to coast
distances_coast_adaptive <- pointDistance(p1 = mesh_adaptive[, c("long", "lat")], 
                                          p2 = coast_longlat, 
                                          lonlat = TRUE)
mesh_adaptive$distance_to_coast <- Rfast::rowMins(distances_coast_adaptive, 
                                                  value = TRUE)
# Visually check the new adaptive mesh
p_depth <- ggplot(data = mesh_adaptive) +
  geom_point(mapping = aes(x = long, y = lat, colour = depth)) + 
  labs(colour = "Depth") + coord_equal()
# + geom_line(data = coast_longlat, mapping = aes(x = long, y = lat, size = "Coastline (ocedata)"), 
#           linetype = 1) + scale_size_manual(values = 1)
p_distance_coast<- ggplot(data = mesh_adaptive) +
  geom_point(mapping = aes(x = long, y = lat, colour = distance_to_coast)) + 
  labs(colour = "D2Coast") + coord_equal()
# + geom_line(data = coast_longlat, mapping = aes(x = long, y = lat, size = "Coastline (ocedata)"), 
#           linetype = 1) + scale_size_manual(values = 1)
p_distance_DASAR<- ggplot(data = mesh_adaptive) +
  geom_point(mapping = aes(x = long, y = lat, colour = to_nearest_DASAR)) + 
  labs(colour = "D2DASAR") + coord_equal()
# + geom_line(data = coast_longlat, mapping = aes(x = long, y = lat, size = "Coastline (ocedata)"), 
#           linetype = 1) + scale_size_manual(values = 1)
p_area <- ggplot(data = mesh_adaptive) +
  geom_point(mapping = aes(x = long, y = lat, colour = area)) + 
  labs(colour = "Area") + coord_equal()
# + geom_line(data = coast_longlat, mapping = aes(x = long, y = lat, size = "Coastline (ocedata)"), 
#           linetype = 1) + scale_size_manual(values = 1)
gridExtra::grid.arrange(p_depth, p_distance_coast, p_distance_DASAR, p_area,
                        ncol = 1, nrow = 4, 
                        top = "The adaptive grid with four different colourings (depth, distance_to_coast, distance_to_DASAR, and area)")

# Write .csv file =========================================================  ##
write.csv(x = mesh_adaptive,
          file = "Data/grid_adaptive_levels=1_maxD2C=100k_maxD2A=100k_area=25966.92_n=2529.csv",
          row.names = FALSE)
