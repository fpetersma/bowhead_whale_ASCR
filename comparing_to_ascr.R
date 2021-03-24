# Load libraries ===============================================================
library(ascr)
library(dplyr)
library(secr)
library(reshape2)
library(readr)

trunc_level <- 96
WITH_NOISE <- FALSE
SINGLE_SL <- TRUE
min_detections <- 1

# Load the files ===============================================================
mesh <- read.csv("Data/grid_adaptive_levels=1_maxD2C=100k_maxD2A=100k_area=21079.44_n=375.csv")
detections <- read.csv("Data/detections_31-08-2010_all.csv")
bearings <- read.csv("Data/bearings_31-08-2010_all.csv")
received_levels <- read.csv("Data/received_levels_31-08-2010_all.csv")
noise_call <- read.csv("Data/noise_31-08-2010_all.csv") # 99% quantile of noise is 95.6
DASAR <- as.data.frame(read_tsv("Data/DASARs.txt"))

DASAR$long <- -abs(DASAR$long) # Turn longitude from westing to easting (more common)
detectors <- DASAR[DASAR$year == 2010 & DASAR$site == 5, c("long", "lat")] 

#===========# Run next section to filter on received levels #==================#

if (WITH_NOISE) {
  sufficient_snr <- received_levels - noise_call >= trunc_level # get indices 
  sufficient_snr[is.na(sufficient_snr)] <- FALSE
  
  detections <- sufficient_snr * 1
  enough_dets <- Rfast::rowsums(detections) >= min_detections
  
  # only keep calls with enough detections, after truncation at 15dB snr
  sufficient_snr <- sufficient_snr[enough_dets, ]
  detections <- detections[enough_dets, ]
  
  bearings <- as.matrix(bearings)[enough_dets, ]
  bearings[!sufficient_snr] <- NA
  
  noise_call <- as.matrix(noise_call)[enough_dets, ]
  
  received_levels <- as.matrix(received_levels)[enough_dets, ]
  received_levels[!sufficient_snr] <- NA
} else if (!WITH_NOISE) {
  sufficient_rl <- received_levels  >= trunc_level # get indices 
  sufficient_rl[is.na(sufficient_rl)] <- FALSE
  
  detections <- sufficient_rl * 1
  enough_dets <- Rfast::rowsums(detections) >= min_detections
  
  # only keep calls with enough detections, after truncation at 15dB snr
  sufficient_rl <- sufficient_rl[enough_dets, ]
  detections <- detections[enough_dets, ]
  
  bearings <- as.matrix(bearings)[enough_dets, ]
  bearings[!sufficient_rl] <- NA
  
  noise_call <- as.matrix(noise_call)[enough_dets, ]
  
  received_levels <- as.matrix(received_levels)[enough_dets, ]
  received_levels[!sufficient_rl] <- NA
}


# Extract information and modify for 'ascr' ====================================
det_hist <- detections

bearings_hist <- bearings 
bearings_rad <- circular::conversion.circular(bearings_hist, units = "radians")
bearings_rad[is.na(bearings_rad)] <- 0

rl_hist <- received_levels
rl_hist[is.na(rl_hist)] <- 0

# snr <- dat_sim$received_levels - dat_sim$noise_call
# snr[is.na(snr)] <- 0

detector_IDs <- DASAR %>% 
  dplyr::filter(year == 2010, site == 5) %>% 
  dplyr::select(ID)
detector_IDs <- unlist(detector_IDs)

# combine the three matrices as one data frame and add id column

capthist <- cbind(id = 1:nrow(det_hist), as.data.frame(det_hist))
colnames(capthist)[2:7] <- detector_IDs
capthist_molten <- melt(capthist, id.vars = "id")

bearhist <- cbind(id = 1:nrow(bearings_rad), as.data.frame(bearings_rad))
colnames(bearhist)[2:7] <- detector_IDs
bearhist_molten <- melt(bearhist, id.vars = "id")

rlhist <- cbind(id = 1:nrow(rl_hist), as.data.frame(rl_hist))
colnames(rlhist)[2:7] <- detector_IDs
rlhist_molten <- melt(rlhist, id.vars = "id")

capt_bear_rl_hist <- cbind(capthist_molten, bear = bearhist_molten$value,
                           rl = rlhist_molten$value)
final_capthist <- data.frame(Session = as.numeric(1), 
                             ID = capt_bear_rl_hist$id,
                             Occasion = as.numeric(1),
                             Detector = capt_bear_rl_hist$variable, 
                             # bear = capt_bear_rl_hist$bear,
                             Signal = capt_bear_rl_hist$rl) #,
#Bearing = capt_bear_snr_hist$bear)
finalCH <- na.omit(final_capthist)


# Create trap object
final_traps <- read.traps(data = data.frame(x = detectors[, 1],
                                            y = detectors[, 2]),
                          detector = "signal")
rownames(final_traps) <- detector_IDs
use <- matrix(1, nr = nrow(final_traps), nc = 1)
# fill it in from your data

usage(final_traps) <- use

correct_capthist <- make.capthist(captures = finalCH,
                                  traps = final_traps, 
                                  fmt = "trapID",
                                  cutval = 96)
ascr_capthist <- list(bincapt = as.matrix(det_hist),
                      ss = as.matrix(rl_hist))


mask <- as.matrix(data.frame(x = mesh$long, 
                             y = mesh$lat))
attr(mask, 'area') <- 56.2
attr(mask, 'buffer') <- 100000
attr(mask, 'dimnames') <- list(1:nrow(mesh), c("x","y"))
mask_secr <- convert.mask(mask = mask)
# convert.mask(mask = example.data$mask)
# traps <- cbind(Detector = detector_IDs, detectors)


a <- fit.ascr(capt = ascr_capthist, 
              traps = final_traps, 
              mask = mask_secr, 
              detfn = "ss", 
              sv = list(D = 1,
                        b0.ss = 150,
                        b1.ss = 1,
                        sigma.ss = 1),
              ss.opts = list(cutoff = 96))















# Create correct object
write.table(final_capthist, file = "capthist_test.txt", 
            row.names = F, 
            col.names = F)
write.table(traps, file = "trap_test.txt", 
            row.names = F, 
            col.names = F)

write.capthist(finalCH, "temp")
correct_capthist <- read.capthist("capthist_test.txt", 
                                  "trap_test.txt", 
                                  detector = "signal", 
                                  cutval = 96)

CH_ascr <- convert.capt.to.ascr(correct_capthist)
traps_ascr <- traps(correct_capthist)

# Fit model


# Test script
write.capthist(signalCH, "temp") ## export data for demo
tempCH <- read.capthist("tempcapt.txt", "temptrap.txt", detector = "signal", cutval = 52.5, )
ovenbird.capt <- convert.capt.to.ascr(signalCH)
ovenbird.traps <- traps(signalCH)
