plotDensity <- function(d) {
  # Takes either a matrix of densities with x and y coordinates, or an
  # object of class "bw_ascr_fit" containing a density matrix, and creates a 
  # density/heat map. 
  
  # Load necessary library 'lattice'
  library(lattice)
  library(raster)
  library(ggplot2)
  
  if (class(d) == "bwASCR_model") {
    d <- d$density
  }

  # ggplot(data = d, mapping = aes(x = long, y = lat, fill = density)) +
  #   geom_tile()
  
  # Create the heat map
  par(mar=c(3,4,2,2))
  levelplot(density ~ long * lat,
            data = d,
            pretty = TRUE,
            xlab = "Longitude easting",
            ylab = "Latitude northing",
            col.regions = heat.colors(100)[length(heat.colors(100)):1],
            main = "")

}
