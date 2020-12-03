####  --------------------------------------------------------------------  ####
####  Exploring the implementation of lat-long instead of UTM coordinates.  ####
####  --------------------------------------------------------------------  ####
####  The problem with this study is that my study are covers multiple UTM  ####
####  zones (zone 6 and 7). These zones are incompatible with each other,   ####
####  and using just one projection for the entire study area will lead to  ####
####  significant biases in angles and/or distances calculated. Instead, I  ####
####  want to work with WGS84 lat-long coordinates. How to do it?           ####  
####                                                                        ####
####  To calculate the angles between points, I will use the function       ####
####  maptools::gzAzimuth().                                                ####
####                                                                        ####
####  To calculate the distances between points, I will use the function    ####
####  geosphere::distm()/raster::pointDistance()                            ####
####                                                                        ####
####  In this script, I will explore how this works and check if it gives   #### 
####  the same distances as the UTM coordinates do.                         ####
####  --------------------------------------------------------------------  ####

#-------------------------# Load the libraries #-------------------------------#
library(tidyverse)
library(raster)
library(geosphere)
library(maptools)
library(marmap)
library(readr)
library(oce)
library(rgdal)
library(circular)
# library(swfscMisc)

#----------------------------# Load the data #---------------------------------#
# Load the detector data for site 5 from 2010
DASAR <- as.data.frame(read_tsv("Data/DASARs.txt"))
DASAR <-  DASAR[DASAR$year == 2010 & DASAR$site == 5, ]
detectors <- DASAR[, c("utmx", "utmy", "long", "lat")]
detectors$long <- -1 * detectors$long # Set longitude to easting

# Get the spatial information from NOAA
lon <- c(-144, -143) 
lat <- c(69.5, 71)

# Extract bathymetry data from NOAA
bathy <- getNOAA.bathy(lon1 = lon[1], lon2 = lon[2], 
                       lat1 = lat[1], lat2 = lat[2],
                       resolution = 4, keep = TRUE)
raster_longlat <- as.raster(bathy)

# Convert to the right UTM coordinates and extract those
raster_utm <- projectRaster(raster_longlat, crs = "+proj=utm +zone=7")

# Use plot() for bathy, raster_longlat and raster_utm -- it will us it the 
# correct plot function, which makes it easy to check if coordinates are correct.
plot(bathy)
plot(raster_longlat)
plot(raster_utm)

# If these are correct, proceed to converting to matrices
# longlat
coord <- coordinates(raster_longlat)
matrix_longlat <- as.matrix(raster_longlat)
rownames(matrix_longlat) <- rev(sort(unique(coord[, 2])))
colnames(matrix_longlat) <- sort(unique(coord[, 1]))

# UTM
coord <- coordinates(raster_utm)
matrix_utm <- as.matrix(raster_utm)
rownames(matrix_utm) <- rev(sort(unique(coord[, 2])))
colnames(matrix_utm) <- sort(unique(coord[, 1]))

# Set all land to NA
matrix_longlat[matrix_longlat > 0] <- NA
matrix_utm[matrix_utm > 0] <- NA

#---------------# Get the coastlines using oce #-------------------------------#
data(coastlineWorldFine, package = "ocedata")
coast <- na.omit(as.data.frame(coastlineWorldFine@data))
coast_longlat <- coast[coast$longitude >= lon[1] &
                        coast$longitude <= lon[2] & 
                        coast$latitude >= lat[1] &
                        coast$latitude <= lat[2], ]
colnames(coast_longlat) <- c("long", "lat")

coast_utm <- project(as.matrix(coast_longlat), 
                     proj = "+proj=utm +zone=7 +datum=WGS84")
colnames(coast_utm) <- c("utmx", "utmy")

#-------------# Calculate the distances and compare #--------------------------#
# Using UTM
melt_utm <- melt(matrix_utm, varnames = c("utmy", "utmx"), 
                 value.name = "depth")
melt_utm <- na.omit(melt_utm)
distances_utm <- pointDistance(p1 = melt_utm[, c("utmx", "utmy")], 
                               p2 = DASAR[, c("utmx", "utmy")], 
                               lonlat = FALSE)

# Using longlat
melt_longlat <- melt(matrix_longlat, varnames = c("lat", "long"), 
                    value.name = "depth")
melt_longlat <- na.omit(melt_longlat)
distances_longlat <- pointDistance(p1 = melt_longlat[, c("long", "lat")], 
                                   p2 = DASAR[, c("long", "lat")], 
                                   lonlat = TRUE)

#----------------# Calculate the angles and compare #--------------------------#
# Using UTM
bearings_utm <- apply(detectors[, c("utmx", "utmy")], 1, function(det) {
  x <- melt_utm$utmx - det[1]
  y <- melt_utm$utmy - det[2]
  return(coord2rad(x, y, control.circular = list(template = "geographics", 
                                                 modulo = "2pi",
                                                 units = "degrees")))
})
hist(bearings_utm)

# Using longlat. bearingRhumb() calculates a rhumb line/loxodrome, which is a path 
# with constant bearing relative to true north. This sounds right, but probably isn't, 
# as a rhumbline is not the shortest distance. However, we do use the shortest 
# distances when calculating distance, so I will use bearing(). 
bearings_longlat <- t(apply(melt_longlat[, c("long", "lat")], 1, function(grid_point) {
  bear <- geosphere::bearing(p1 = as.matrix(detectors[, c("long", "lat")]),
                             p2 = grid_point)
}))

bearings_longlat <- as.circular(bearings_longlat, 
                                control.circular = list(template = "geographics", 
                                                        modulo = "2pi", 
                                                        units = "degrees"))

# Evaluate the similarity in distribution visually.
ggplot(data = rbind(data.frame(bearing = as.numeric(bearings_longlat), 
                               unit = "long-lat"),
                    data.frame(bearing = as.numeric(bearings_utm),
                               unit = "utm"))) +
  geom_histogram(mapping = aes(x = bearing, fill = unit), position = "identity",
                 alpha = 0.4)
