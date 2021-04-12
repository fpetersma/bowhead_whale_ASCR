
library(circular)
library(tidyverse)

detectors <- as.data.frame(readr::read_tsv("Data/DASARs.txt")) %>% 
  filter(year == 2010 & site == 5) %>% 
  select(utmx, utmy, lat, long, ID) 

bearings <- read.csv("~/University of St Andrews/PhD/Bowhead Whales/JABES paper/Data/bearings_31-08-2010_all.csv",
                     na.strings = -1, col.names = detectors$ID)


bearings_deg <- circular(bearings, 
                         units = "degrees",
                         template = "geographics",
                         # rotation = "counter",
                         # zero = 0.5*pi,
                         ) # 'geographics' means clockwise rotation with 0 at north


plot(bearings_deg[1, 5])
plot(bearings_deg[1, 6])

bearings_rad <- conversion.circular(bearings_deg)

plot(bearings_rad[1, 5])
plot(bearings_rad[1, 6])

# NEW:
plotBearings <- function(link) {
  pos_det <- !is.na(link)
  pos_detectors <- detectors[pos_det, ]
  x <- -pos_detectors$long
  y <- pos_detectors$lat
  rad <- link[pos_det]
  
  # Create map
  plot(-detectors$long, detectors$lat, pch = 21, bg = "white",
       col = "red", cex = 1, xlim = c(-144, -142.5),
       ylim = c(70, 70.7))
  points(x, y, pch = 21, bg = "red",
         col = "red", cex = 1)
  
  # Add arrows
  arrows(x, y, x + sin(rad) * 0.1, y + cos(rad) * 0.1, angle = 15)
}

plotBearings(link = bearings_rad[1, ])

# Plot several
index <- sample(1:nrow(bearings_rad), size = 100, replace = FALSE)
selection <- bearings_rad[index, ]


for (i in 1:nrow(selection)) {
  readline(prompt="Press [enter] to create next plot")
  linkage <- selection[i, ]
  plotBearings(linkage)
  title(main = paste("Row ID: ",index[i]))
}


