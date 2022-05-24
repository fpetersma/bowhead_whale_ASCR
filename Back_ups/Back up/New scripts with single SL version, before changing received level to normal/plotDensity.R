plotDensity <- function(d) {
  # Takes either a matrix of densities with x and y coordinates, or an
  # object of class "bw_ascr_fit" containing a density matrix, and creates a 
  # density/heat map. 
  
  # Load necessary library 'lattice'
  library("lattice")
  
  if (class(d) == "bw_ascr_fit") {
    d <- d$density
  }
  
  # Create the heat map
  # par(mar=c(3,4,2,2))
  levelplot(density ~ x * y, 
            data = d, 
            pretty = TRUE,
            xlab = "UTM easting", 
            ylab = "UTM northing",
            col.regions = heat.colors(100)[length(heat.colors(100)):1],
            main = "")
}