# Load the data
load(file = "Data/JABES paper/100_simulated_datasets_N=4471.RData")
dat <- bw_list[[1]]
attach(dat)

source_levels <- 120:180
b <- 100
random_noise_single <- noise_random[sample(1:1000, size = b), 1] 

df <- function(sl, d, noise) {
  SNR <- sl - 15 * log(d, base = 10) - noise 
  p <- 0.8 * pnorm(SNR, mean = 15, sd = 7) * dnorm(sl, mean = 154, sd = 9)
  p
}

prob <- mean(sapply(random_noise_single, function(noise) {
  sum(df(source_levels, 150000, noise))
}))

p <- 0.8 * pnorm(150 - 15 * log(200000, base = 15) - 70, mean = 20, sd = 7)



### HO, this is for a single detector, but that is not what we are interest in. 
### We are interested in having at least TWO detections.
source(file = "Scripts/Bowhead Whales/hidden_functions.R")

distances <- matrix(c(150000, 145000, 140000, 142000, 149000, 138000),
                    nrow = 1)
distances <- c(200000, 210000, 215000, 211000, 218000, 280000)

probs <- sapply(distances, function(d) {
  prob <- mean(sapply(random_noise_single, function(noise) {
    return(sum(df(100:200, d, noise)))
  }))
  return(prob)
})
.detected(probs = matrix(probs, nrow = 1), min_no_detections = 2)


### what if we turn this around, because I think I might have violated some mathematical rules here...
source_levels <- 120:180
random_noise_single <- noise_random[sample(1:1000, size = 100), 1]
distances <- rep(300000, 6)
distances <- c(200000, 210000, 215000, 211000, 218000, 208000)

expected_p. <- mean({
  Rfast::colsums(apply(noise_random, 1, function(c) {
    sl_m <- matrix(rep(source_levels, each = 6), byrow = TRUE, ncol = 6)
    c_m <- matrix(c, nrow = length(source_levels), ncol = 6, byrow = TRUE)
    distances_m <- matrix(distances, nrow = length(source_levels), ncol = 6, byrow = TRUE)
    
    E_snr <- sl_m - 15 * log(distances_m, base = 10) - c_m
    p <- 0.7 * pnorm(E_snr, mean = 20, sd = 5) 
    
    p. <- .detected(p, 2) * dnorm(sl_m[, 1], mean = 140, sd = 5)
  }))
})

# Create a map
# Load the DASAR locations
year <- 2010
site <- 5

DASAR <- as.data.frame(readr::read_tsv("Data/DASARs.txt"))
DASAR <-  DASAR[DASAR$year == year & DASAR$site == site, ]
detectors <- subset(DASAR, select = c("long", "lat"))

mesh <- read.csv("Data/JABES paper/grid_adaptive_levels=3_bounds=10k-60k_maxD2C=100k_maxD2A=200k_area=44145.6_n=180+139+146=465.csv")
########################### Create mask ########################################
cov_density <- mesh
cov_density$depth <- abs(cov_density$depth) # make depth positive 
cov_density$depth[cov_density$depth == 0] <- 0.01
cov_density[["depth2"]] <-  cov_density$depth ^ 2
cov_density[["logdepth"]] <- log(cov_density$depth)
cov_density[["logdepth2"]] <- cov_density$logdepth ^ 2
cov_density[["distance_to_coast2"]] <- cov_density$distance_to_coast ^ 2
# Create a matrix of distances
distances <- pointDistance(p1 = cov_density[, c("long", "lat")], 
                           p2 = detectors[, c("long", "lat")], 
                           lonlat = TRUE)

prob_grid <- apply(distances, 1, function(d) {
  expected_p. <- mean({
    Rfast::colsums(apply(noise_random[1:100, ], 1, function(c) {
      sl_m <- matrix(rep(source_levels, each = 6), byrow = TRUE, ncol = 6)
      c_m <- matrix(c, nrow = length(source_levels), ncol = 6, byrow = TRUE)
      d_m <- matrix(d, nrow = length(source_levels), ncol = 6, byrow = TRUE)
      
      E_snr <- sl_m - 15 * log(d_m, base = 10) - c_m
      p <- 0.7 * pnorm(E_snr, mean = 20, sd = 5) 
      
      p. <- .detected(p, 2) * dnorm(sl_m[, 1], mean = 154, sd = 9)
    }))
  })
})

new <- cbind(grid[, c("long", "lat")], probs = prob_grid)
new <- reshape2::dcast(data = new, lat ~ long, value.var = "probs")
rownames(new) <- new[, 1]
new <- new[, -1]
plotly::plot_ly(z = ~ as.matrix(new)) %>% plotly::add_surface()
