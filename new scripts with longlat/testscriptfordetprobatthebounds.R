## If we assume a 

source_levels <- 100:200
b <- 100
random_noise_single <- noise_random[sample(1:1000, size = b), 1] 

df <- function(sl, d, noise) {
  SNR <- sl - 15 * log(d, base = 10) - noise 
  p <- 0.8 * pnorm(SNR, mean = 15, sd = 7) * dnorm(sl, mean = 154, sd = 8)
  p
}

prob <- mean(sapply(random_noise_single, function(noise) {
  sum(df(100:200, 150000, noise))
}))

p <- 0.8 * pnorm(150 - 15 * log(100000, base = 10) - 70, mean = 20, sd = 7)



### HO, this is for a single detector, but that is not what we are interest in. 
### We are interested in having at least TWO detections.
source(file = "Scripts/Bowhead Whales/hidden_functions.R")

distances <- matrix(c(150000, 145000, 140000, 142000, 149000, 138000),
                    nrow = 1)
distances <- c(150000, 145000, 140000, 142000, 149000, 138000)

probs <- sapply(distances, function(d) {
  prob <- mean(sapply(random_noise_single, function(noise) {
    return(sum(df(100:200, d, noise)))
  }))
  return(prob)
})
.detected(probs = matrix(probs, nrow = 1), min_no_detections = 2)


### what if we turn this around, because I think I might have violated some mathematical rules here...
source_levels <- 100:200
random_noise_single <- noise_random[sample(1:1000, size = 100), 1]
distances <- rep(150000, 6)
distances <- c(150000, 145000, 140000, 142000, 149000, 138000)

expected_p. <- mean({
  Rfast::colsums(apply(noise_random[1:100, ], 1, function(c) {
    sl_m <- matrix(rep(source_levels, each = 6), byrow = TRUE, ncol = 6)
    c_m <- matrix(c, nrow = length(source_levels), ncol = 6, byrow = TRUE)
    distances_m <- matrix(distances, nrow = length(source_levels), ncol = 6, byrow = TRUE)
    
    E_snr <- sl_m - 15 * log(distances_m, base = 10) - c_m
    p <- 0.6 * pnorm(E_snr, mean = 20, sd = 10) 
    
    p. <- .detected(p, 2) * dnorm(sl_m[, 1], mean = 154, sd = 8)
  }))
})

# Create a map
# Load the DASAR locations
year <- 2010
site <- 5

DASAR <- as.data.frame(readr::read_tsv("Data/DASARs.txt"))
DASAR <-  DASAR[DASAR$year == year & DASAR$site == site, ]
detectors <- subset(DASAR, select = c("utmx", "utmy"))

grid <- read.csv("Data/Masks/[PT] mask_buffer=60km_mesh=coarse_A=21597576m2_deltaX=4466_deltaY=4836_grid=393.csv")

# Create a matrix of distances
distances <- apply(detectors, 1, function(z) {
  y <- t(grid[, c("utmx", "utmy")]) - z
  y <- y ^ 2
  return(sqrt(Rfast::colsums(y)))
})

prob_grid <- apply(distances, 1, function(d) {
  expected_p. <- mean({
    Rfast::colsums(apply(noise_random[1:100, ], 1, function(c) {
      sl_m <- matrix(rep(source_levels, each = 6), byrow = TRUE, ncol = 6)
      c_m <- matrix(c, nrow = length(source_levels), ncol = 6, byrow = TRUE)
      d_m <- matrix(d, nrow = length(source_levels), ncol = 6, byrow = TRUE)
      
      E_snr <- sl_m - 15 * log(d_m, base = 10) - c_m
      p <- 0.6 * pnorm(E_snr, mean = 20, sd = 10) 
      
      p. <- .detected(p, 2) * dnorm(sl_m[, 1], mean = 154, sd = 8)
    }))
  })
})

new <- cbind(grid[, c("utmx", "utmy")], probs = prob_grid)
new <- reshape2::dcast(data = new, utmy ~ utmx, value.var = "probs")
rownames(new) <- new[, 1]
new <- new[, -1]
plot_ly(z = ~ as.matrix(new)) %>% add_surface()
