####  --------------------------------------------------------------------  ####
####            CREATE GRID WITH DISTANCE-TO-COAST AND DEPTH                ####
####  --------------------------------------------------------------------  ####
####  In this script, I will extract bathemetry data from NOAA and derive   ####
####  the distance to coast using the dist2Line() function.                 ####
####                                                                        ####
####  --------------------------------------------------------------------  ####

#---------------------------# Load libraries #---------------------------------#
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

#--------------------------# Create the map #----------------------------------#

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
                       resolution = 1, keep = TRUE)
raster_longlat <- as.raster(bathy)

# Check area
a <- area(raster_longlat)
a[, 1]
sum(a[]) # should be 204387.3 km2

# If these are correct, proceed to converting to matrices
coord <- coordinates(raster_longlat)
matrix_longlat <- as.matrix(raster_longlat)
rownames(matrix_longlat) <- rev(sort(unique(coord[, 2])))
colnames(matrix_longlat) <- sort(unique(coord[, 1]))

# Set all land to NA
matrix_longlat[matrix_longlat >= 0] <- NA

#------------------# Get the coastlines using oce #----------------------------#
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

# # Below is old and not trustworthy
# distance_to_coast <- dist2Line(p = as.matrix(select(melt_longlat, long, lat)), 
#                                line = as.matrix(coast_longlat))

melt_longlat <- cbind(melt_longlat, distance_to_coast = distance_to_coast)

# To visually check, run ggplot below
ggplot(data = melt_longlat) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  # geom_point(mapping = aes(colour = to_nearest_DASAR, fill = to_nearest_DASAR), show.legend = TRUE) +
  geom_contour_filled(mapping = aes(x = long, y = lat, z = depth), binwidth = 400) +
  geom_line(data = coast_longlat, mapping = aes(x = long, y = lat, size = "Coastline (ocedata)"), 
            linetype = 1) +
  geom_contour(mapping = aes(x = long, y = lat, z = to_nearest_DASAR, colour = "From array\n(per 50km)"), 
               linetype = 2, binwidth = 50000) +
  geom_contour(mapping = aes(x = long, y = lat, z = distance_to_coast, colour = "From coast\n(per 25km)"), 
               linetype = 2, binwidth = 25000) +
  geom_point(data = detectors, mapping = aes(x = long, y = lat, shape = "DASAR array"), 
             colour = "red", size = 2) +
  scale_colour_manual(values = c("red", "black", "red")) +
  scale_size_manual(values = 1) +
  scale_shape_manual(values = 17) +
  labs(fill = "Depth (m)", colour = "Distance contours", size = "", shape = "", 
       y = "Latitude (Northing)", x = "Longitude (Easting)", 
       title = "Map of study area", 
       subtitle = "with depth, distance to coast and distance to DASAR array superimposed") 

#----------------------# Create adaptive mesh #--------------------------------#
buffer_50 <- subset(melt_longlat, subset = to_nearest_DASAR <= 50000)
buffer_100 <- subset(melt_longlat, subset = to_nearest_DASAR <= 100000)
buffer_200_100 <- subset(melt_longlat, subset = to_nearest_DASAR <= 200000 &
                           distance_to_coast <= 100000)

# # Let's try inlabru
# mesh <- INLA::inla.mesh.2d(loc = buffer_200_100[, c(2, 1)], max.edge = 0.5, 
#                            cutoff = 0.5)
# plot(mesh) # Okay, no idea if this is useful, but it's fun =D

#----# Let's try my own method #----#
# Create 200 grid points in first circle, 100 in second circle, and 100 in third circle.
# I created a buffer using three different resolutions, using a coarser resolution
# the farther away fromt the detector array we get. This meant re-running this whole
# script and every time storing only the relevant part of the mesh at the preferred 
# resolution (see circle_1, circle_2, circle_3).
mesh <- buffer_200_100

mesh_raster <- rasterFromXYZ(mesh[, 1:3], 
                             crs = "+proj=longlat +datum=WGS84 +no_defs",
                             res = c(1/60, 1/60)) # res depends on the resolution
# when creating a raster using rasterFRomXYZ(), it adds raster cells that contain 
# NA values for depth. We are not interested in those, so when we check the total 
# area, we only want to sum the values that actually have a depth measurement 
# associated with them.
sum(area(mesh_raster)[][!is.na(values(mesh_raster))]) # values(mesh_raster) gives the depth values

# circle 1 = <=50k, 50k < circle 2 <=100k, circle 3 > 100k
mesh$circle <- ifelse(test = mesh$to_nearest_DASAR <= 50000,
                      yes = 1, no = 2)
mesh$circle[mesh$to_nearest_DASAR > 100000] <- 3

sum(mesh$circle == 1) # 7223 @ res=1 c(1/60, 1/60), 201 @ res=6, c(1/10, 1/10)
sum(mesh$circle == 2) # 11393 @ res=1 c(1/60, 1/60), 120 @ res=10, c(1/6, 1/6)
sum(mesh$circle == 3) # 18935 @ res=1 c(1/60, 1/60), 94 @ res=14, c(7/30, 7/30)

# # No need to run next 3 lines again, as data was saved as .RData.
# circle_1 <- mesh[mesh$circle == 1, ] 
# circle_2 <- mesh[mesh$circle == 2, ]
# circle_3 <- mesh[mesh$circle == 3, ]

# Transform the circles to rasters to get the area for every cell. As the earth
# is not flat (duh), the area per cell varies not just between circles, but also
# between different latitude (and longitudes maybe, but I don't think so). Either 
# way, we need an individual area measurement for every grid cells we integrate 
# over. 
raster_1 <- rasterFromXYZ(circle_1[, c(2, 1, 3)], 
              crs = "+proj=longlat +datum=WGS84 +no_defs",
              res = c(0.1, 0.1))
raster_2 <- rasterFromXYZ(circle_2[, c(2, 1, 3)], 
              crs = "+proj=longlat +datum=WGS84 +no_defs",
              res = c(1/6, 1/6))
raster_3 <- rasterFromXYZ(circle_3[, c(2, 1, 3)], 
              crs = "+proj=longlat +datum=WGS84 +no_defs",
              res = c(7/30, 7/30))

# Here we recreate the melted data we had before. 
df_1 <- na.omit(cbind(coordinates(raster_1), 
                      depth = values(raster_1), 
                      area = area(raster_1)[]))
df_2 <- na.omit(cbind(coordinates(raster_2), 
                      depth = values(raster_2), 
                      area = area(raster_2)[]))
df_3 <- na.omit(cbind(coordinates(raster_3), 
                      depth = values(raster_3), 
                      area = area(raster_3)[]))

# Put it all together
mesh_adaptive <- as.data.frame(rbind(df_1, df_2, df_3))
colnames(mesh_adaptive) <- c("long", "lat", "depth", "area")
sum(mesh_adaptive$area) # 43384.24, which is a little off the true, which is 
                        # 43176.2, which is probably due to a bit of overlap. 
                        # I don't think this is a major issue. 

# mesh_adaptive <- rbind(circle_1, circle_2, circle_3)

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
  labs(colour = "Depth")
  # + geom_line(data = coast_longlat, mapping = aes(x = long, y = lat, size = "Coastline (ocedata)"), 
  #           linetype = 1) + scale_size_manual(values = 1)
p_distance_coast<- ggplot(data = mesh_adaptive) +
  geom_point(mapping = aes(x = long, y = lat, colour = distance_to_coast)) + 
  labs(colour = "D2Coast")
  # + geom_line(data = coast_longlat, mapping = aes(x = long, y = lat, size = "Coastline (ocedata)"), 
  #           linetype = 1) + scale_size_manual(values = 1)
p_distance_DASAR<- ggplot(data = mesh_adaptive) +
  geom_point(mapping = aes(x = long, y = lat, colour = to_nearest_DASAR)) + 
  labs(colour = "D2DASAR")
  # + geom_line(data = coast_longlat, mapping = aes(x = long, y = lat, size = "Coastline (ocedata)"), 
  #           linetype = 1) + scale_size_manual(values = 1)
p_area <- ggplot(data = mesh_adaptive) +
  geom_point(mapping = aes(x = long, y = lat, colour = area)) + 
  labs(colour = "Area")
  # + geom_line(data = coast_longlat, mapping = aes(x = long, y = lat, size = "Coastline (ocedata)"), 
  #           linetype = 1) + scale_size_manual(values = 1)
gridExtra::grid.arrange(p_depth, p_distance_coast, p_distance_DASAR, p_area,
                        ncol = 1, nrow = 4, 
                        top = "The adaptive grid with four different colourings (depth, distance_to_coast, distance_to_DASAR, and area)")

#-----------------------# Write .csv file #------------------------------------#
write.csv(x = mesh_adaptive, 
          file = "Data/JABES paper/grid_adaptive_maxD2C=100k_maxD2A=200k.csv", 
          row.names = FALSE)
