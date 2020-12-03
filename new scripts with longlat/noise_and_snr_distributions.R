################################################################################
####  "noise_and_snr_distributions.R"                                       ####
####  --------------------------------------------------------------------  ####
####  Explore the distribution of noise and SNR. The main challenge is to   ####
####  to find out what the distribution of noise levels is, and whether     ####
####  differ whether a call was detected or not. Also, I need to see how I  ####
####  handle negative SNR values.                                           ####
####  --------------------------------------------------------------------  ####
####  By: Felix Petersma                                                    ####
####  When: 25/08/2020                                                      ####
################################################################################

#### Load libraries ############################################################

library(tidyverse)

#### Read in data ##############################################################

detections <- read_csv("../JABES paper/Data/detections_31-08-2020_successful-loc.csv")
noise <- read_csv("../JABES paper/Data/noise_31-08-2020_successful-loc.csv")
snr <- read_csv("../JABES paper/Data/snr_31-08-2020_successful-loc.csv")

detections <- as.matrix(detections)
noise <- as.matrix(noise)
snr <- as.matrix(snr)

noise_flat <- data.frame(noise = as.vector(noise), detected = as.factor(as.vector(detections)))

ggplot(noise_flat) + geom_density(aes(x = noise, group = detected, fill = detected), alpha = 0.4) + 
  geom_density(aes(x = noise), alpha = 0.7, fill = "black") + theme_bw() + 
  ggtitle("Distribution of all noise, noise for detections (blue) and noise for non-detections (red), for site 5, 2010.")

# Noise seems to be lower for calls where the sound was not detected... How could that be?

hist(noise)
hist(noise[as.logical(detections)])
hist(noise[!as.logical(detections)])
