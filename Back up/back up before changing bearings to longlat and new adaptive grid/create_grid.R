################################################################################
####            CREATE GRID WITH DISTANCE-TO-COAST AND DEPTH                ####
####  --------------------------------------------------------------------  ####
####  In this script, I will extract bathemetry data from NOAA and derive   ####
####  the distance to coast using the dist2Line() function.                 ####
####                                                                        ####
################################################################################

#### LOAD LIBRARIES ############################################################

library(oceanmap)
library(readr)
library(reshape2)
library(geosphere)
library(rgdal)

#### CREATE THE MAP ############################################################

# Load the DASAR locations
year <- 2010
site <- 5

DASAR <- as.data.frame(read_tsv("Data/DASARs.txt"))
DASAR <-  DASAR[DASAR$year == year & DASAR$site == site, ]
detectors <- DASAR[, c("utmx", "utmy")]

x_rangeDASAR <- range(detectors$utmx)
y_rangeDASAR <- range(detectors$utmy)

# Site 5 is in zone 7, site 3 is in zone 6
lon <- c(-144, -138) # The entire width of UTM zone 7
lat <- c(69.5, 72) # Hopefully enough

# Extract bathymetry data from NOAA
bathy <- get.bathy(lon=lon, lat=lat, keep = TRUE, # Is not working atm, so load earlier extraction
                   visualize = FALSE, grid = FALSE, resolution = 0.1)
load(file = "bathy_lon-144--138.lat69.5-71_res.0.1min.dat")
bathy <- h
rm(h)
# Convert to the correct UTM coordinates and extract those
bathy_utm <- projectRaster(bathy, crs = "+proj=utm +zone=7 +ellps=WGS84")
coord <- coordinates(bathy_utm)

# Extract bathymetry data and coordinates as matrix
bathy_matrix <- as.matrix(bathy_utm)
rownames(bathy_matrix) <- round(rev(sort(unique(coord[, 2]))))
colnames(bathy_matrix) <- round(sort(unique(coord[, 1])))

# Melt the data
molten_bathy <- melt(bathy_matrix, varnames = c("utmy", "utmx"), 
                     value.name = "depth")
molten_bathy <- molten_bathy[, c("utmx", "utmy", "depth")]

# remove all land values (= NaN values)
molten_bathy <- na.omit(molten_bathy)
molten_bathy[["dist_to_coast"]] <- NA

# Find the coastline
coastline <- t(sapply(unique(molten_bathy$utmx), function(x) {
  y <- molten_bathy[molten_bathy$utmx == x, "utmy"]
  return(c("utmx" = x, "utmy" = min(y, na.rm = TRUE) - 1))
}))

# save coast line as .csv
# write_csv(as.data.frame(coastline), "Data/coastline.csv") # Already done

# Add the distance to the nearest DASAR
system.time({
  molten_bathy <- cbind(molten_bathy, "distance_to_DASARs" = apply(molten_bathy, 1, function(point) {
    diff <- as.matrix(detectors[, c("utmx", "utmy")]) -
      matrix(as.numeric(point[c("utmx", "utmy")]), nrow = nrow(detectors), ncol = 2, byrow = TRUE)
    
    distances <-  apply(diff, 1, function(d) {
      sqrt(sum(d ^ 2))
    })
    return(min(distances))
  }))
})

temp <- dcast(molten_bathy, formula = utmy ~ utmx, value.var = "distance_to_DASARs")
rownames(temp) <- temp$utmy
temp <- temp[, -1]
plotly::plot_ly(y = as.numeric(rownames(temp[, 1:500])), 
                x = as.numeric(colnames(temp[, 1:500])),
                z = as.matrix(temp[,1:500])) %>% 
  plotly::add_surface()

image(x = molten_bathy$utmx[1:10], 
      y = molten_bathy$utmy[1:10], 
      z = molten_bathy$distance_to_DASARs[1:100])

# Plot the distances to check
library(ggplot2)
with(molten_bathy, plot(x = utmx, y = utmy, fill = distance_to_DASARs))
ggplot()

# Find distance to the coastline
molten_bathy_longlat <- project(as.matrix(molten_bathy[, c("utmx", "utmy")]), 
                                proj = "+proj=utm +zone=7, +ellps=WGS84", 
                                inv = TRUE)
colnames(molten_bathy_longlat) <- c("long", "lat") # bathy in longlat
coastline_longlat <- project(as.matrix(coastline[, c("utmx", "utmy")]), 
                             proj = "+proj=utm +zone=7, +ellps=WGS84", 
                             inv = TRUE)
colnames(coastline_longlat) <- c("long", "lat") # coastline in longlat
system.time({
  dist <- dist2Line(p = molten_bathy_longlat, line = coastline_longlat)
})


#### Save maps as .csv files ###################################################