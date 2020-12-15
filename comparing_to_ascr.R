# Load libraries ===============================================================
library(ascr)
library(dplyr)
library(secr)
library(reshape2)

# Load the files ===============================================================
load("C:/Users/felix/OneDrive - University of St Andrews/Documents/University of St Andrews/PhD/Bowhead Whales/temp_data_object.RData")

# Extract information and modify for 'ascr' ====================================
det_hist <- dat_sim$det_hist

bearings <- dat_sim$bearings
bearings_rad <- circular::conversion.circular(bearings, units = "radians")
bearings_rad[is.na(bearings_rad)] <- 0

snr <- dat_sim$received_levels - dat_sim$noise_call
snr[is.na(snr)] <- 0

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

snrhist <- cbind(id = 1:nrow(snr), as.data.frame(snr))
colnames(rlhist)[2:7] <- detector_IDs
snrhist_molten <- melt(snrhist, id.vars = "id")

capt_bear_snr_hist <- cbind(capthist_molten, bear = bearhist_molten$value,
                            snr = snrhist_molten$value)
final_capthist <- data.frame(Session = 1, 
                             ID = capt_bear_rl_hist$id,
                             Occasion = 1,
                             Detector = capt_bear_snr_hist$variable, 
                             Signal = capt_bear_snr_hist$snr) #,
                             #Bearing = capt_bear_snr_hist$bear)


traps <- cbind(Detector = detector_IDs, detectors)

# Create correct object
write.table(final_capthist, file = "capthist_test.txt", row.names = F, col.names = F)
write.table(traps, file = "trap_test.txt", row.names = F, col.names = F)

correct_capthist <- read.capthist("capthist_test.txt", "trap_test.txt", detector = "signal", cutval = 15)

CH_ascr <- convert.capt.to.ascr(correct_capthist)
traps_ascr <- traps(correct_capthist)



# Test script
write.capthist(signalCH, "temp") ## export data for demo
tempCH <- read.capthist("tempcapt.txt", "temptrap.txt", detector = "signal", cutval = 52.5, )
ovenbird.capt <- convert.capt.to.ascr(signalCH)
ovenbird.traps <- traps(signalCH)
