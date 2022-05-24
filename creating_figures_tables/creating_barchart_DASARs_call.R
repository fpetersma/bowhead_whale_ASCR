##### Script to load the raw detection data from Greeneridge Inc.          #####
##### We focus on the data from Site  5 from 2010.                    #####

library(readr)
library(dplyr)

########################## Load the .tsv files #################################
folder_loc <- "../Export all Aaron/Export_Automated_2010_Initial/2010 [no site 2]/" 
files <- list.files(folder_loc)

site <- 5 # Change this for different sites
n_detectors <- 6

index <- substr(files, 8, 8) == site

files <- files[index] # Only keep files belonging to specific 

# files <- files[49:64] # ONLY USE THIS TO GET SPECIFIC DAYS

colnames_1_25 <- c("time_(AKDT)",
                   "datetime_(AKDT)",
                   "time_(UTC)",
                   "datetime_(UTC)",
                   "UTM_easting",
                   "UTM_northing",
                   "UTM_zone",
                   "latitude_degrees",
                   "longitude_degrees",
                   "distance",
                   "call_type",
                   "area_ellipse_90",
                   "major_axis_(m)",
                   "minor_axis_(m)",
                   "bearing_semi_major",
                   "angle_semi_major_(radians)",
                   "nr_DASARs_used", 
                   "cov_matrix_call_location1",
                   "cov_matrix_call_location2",
                   "cov_matrix_call_location3",
                   "cov_matrix_call_location4",
                   "succesful_location", 
                   "operator_comment", 
                   "operator_ID", 
                   "processing_datetime_(UTC)")

colnames_good <- c("time_(UTC)",
                   "datetime_(UTC)",
                   
                   "UTM_easting",
                   "UTM_northing",
                   "UTM_zone",
                   
                   "distance?",
                   
                   "area_ellipse_90?",
                   "major_axis_(m)?",
                   
                   "bearing_semi_major?",
                   "angle_semi_major_(radians)?",
                   "nr_DASARs_used", 
                   "cov_matrix_call_location1",
                   "cov_matrix_call_location2",
                   "cov_matrix_call_location3",
                   "cov_matrix_call_location4",
                   "succesful_location", 
                   "operator_comment", 
                   "operator_ID", 
                   "processing_datetime_(UTC)")

# Read tsv files and put them in a list
detections_list <- list()

# Before running this loop think about i) how many detectors there are at the site, 
for(i in seq_along(files)) {
  # closer look showed max n_col is 152, so fix at 152.
  n_col <- 19 + 12 * n_detectors # Increasing more
  
  d <- read_tsv(paste0(folder_loc,
                       files[i]), 
                col_names = paste0("V", seq_len(n_col)))
  d <- d[, -ncol(d)] ## again very hacky, but apparently read_tsv adds an empty column... I wonder if everything is broken due to an update in read_tsv...?
  colnames(d)[1:19] <- colnames_good
  cols_detection <- c("DASAR_ID",
                      "location_weight",
                      "brel",
                      "bgrid",
                      "arrival_time_(UTC)",
                      "rms_db_re_1mPa",
                      "SNR_db",
                      "lowest_freq_call",
                      "highest_freq_call",
                      "length_call_(s)",
                      "SE_call_bearings",
                      "kappa_(Von_Mises)")
  cols_leftover <- paste0(cols_detection, "_", rep(1:((n_col-19)/12), each = 12))
  d <- cbind(d, matrix(NA, ncol = n_col - ncol(d), nrow = nrow(d))) # very hacky making sure there are enough columns...
  colnames(d)[20:(n_col)] <- cols_leftover
  detections_list[[i]] <- d #[!is.nan(d$`time_(UTC)`), ] # EXCLUDING ALL NaNs at time_UTC WAS WRoNG!!! BAD DECISION, BUT FIXED NOW
}

# Turn list into on big data.frame
detections <- do.call("rbind", detections_list)

################### Create capture history #####################################

detectors <- as.data.frame(read_tsv("Data/DASARs.txt"))
detectors_2010 <- detectors[detectors$year == 2010 & detectors$site == site, ]

W <- apply(detections, 1, function(w) {
  return(as.numeric(detectors_2010$ID %in% w))
})
W <- data.frame(t(W), stringsAsFactors = FALSE)
colnames(W) <- detectors_2010$ID

head(W)

dets_involved <- as.data.frame(table(rowsums(as.matrix(W))))

library(Rfast)
library(ggplot2)


## load the other data
detections_multi <- read.csv("Data/all-data_all.csv")

n_detections <- c(table(rowsums(as.matrix(W)))[1], table(detections_multi$n_sensors))
dets_involved <- as.data.frame(n_detections)
dets_involved$percentages <- round(dets_involved$n_detections / 
                                     sum(dets_involved$n_detections), 3) * 100
dets_involved$n_DASAR <- 1:6

# hist(rowsums(as.matrix(W)))
ggplot(data = dets_involved, mapping =  aes(x = n_DASAR, y = percentages / 100)) + 
  theme_classic() +
  geom_bar(stat = "identity", fill = "grey75") +
  labs(x = "# DASARs involved", y = "Proportions of detected calls") +
  geom_text(aes(label = paste0(percentages, "%")), 
            position=position_dodge(width = 0.9), vjust = -0.25) +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())



## Everything belwo here is not relevant, 
## but was copied from [PT] Raw data editing.R.



# 
# 
# 
# 
# 
# 
# 
# 
# ################## Select relevant information in dets #########################
# 
# selection <- detections[, paste0(c("DASAR_ID_",
#                                    "bgrid_", 
#                                    "arrival_time_(UTC)_", 
#                                    "rms_db_re_1mPa_",
#                                    "SE_call_bearings_",
#                                    "kappa_(Von_Mises)_"), rep(1:n_detectors, each = 6))]
# det_bearings <- detections[, paste0(c("DASAR_ID_",
#                                       "location_weight_",
#                                       "brel_",
#                                       "bgrid_", 
#                                       "SE_call_bearings_",
#                                       "kappa_(Von_Mises)_"), rep(1:n_detectors, each = 6))]
# only_IDs <- as.matrix(detections[, paste0("DASAR_ID_", rep(1:n_detectors, each = 1))])
# only_bearings <- as.matrix(detections[, paste0("bgrid_", rep(1:n_detectors, each = 1))])
# only_weights <- as.matrix(detections[, paste0("location_weight_", rep(1:n_detectors, each = 1))])
# 
# # bearings_weight_1 <- matrix(NA, nrow = nrow(det_bearings), ncol = nrow(detectors_2010))
# # colnames(bearings_weight_1) <- detectors_2010$ID
# # for (detector in detectors_2010$ID) {
# #   index <- which(only_IDs == detector, arr.ind = TRUE)
# #   valid_bearings_index <- index[as.numeric(only_weights[index]) > 0.99, ] # Only save indices for weight > 0.99
# #   bearings_weight_1[valid_bearings_index[, 1], detector] <- only_bearings[valid_bearings_index]
# # }
# # rm(index, valid_bearings_index)
# 
# bearings_all <- matrix(NA, nrow = nrow(det_bearings), ncol = nrow(detectors_2010))
# colnames(bearings_all) <- detectors_2010$ID
# for (detector in detectors_2010$ID) {
#   index <- which(only_IDs == detector, arr.ind = TRUE)
#   bearings_all[index[, 1], detector] <- only_bearings[index]
# }
# rm(index)
# 
# # Set bearings with weight < 1 to NA
# bearings_successful <- bearings_all
# bearings_successful[detections$succesful_location != "successful", ] <- NA
# 
# weights <- matrix(NA, nrow = nrow(det_bearings), ncol = nrow(detectors_2010))
# colnames(weights) <- detectors_2010$ID
# for (detector in detectors_2010$ID) {
#   index <- which(only_IDs == detector, arr.ind = TRUE)
#   weights[index[, 1], detector] <- only_weights[index]
# }
# rm(index)
# 
# brels <- matrix(NA, nrow = nrow(det_bearings), ncol = nrow(detectors_2010))
# colnames(brels) <- detectors_2010$ID
# for (detector in detectors_2010$ID) {
#   index <- which(det_bearings == detector, arr.ind = TRUE)
#   index[, 2] <- index[, 2] + 2
#   brels[index[, 1], detector] <- as.matrix(det_bearings)[index]
# }
# rm(index)
# 
# von_mises <- matrix(NA, nrow = nrow(det_bearings), ncol = nrow(detectors_2010))
# colnames(von_mises) <- detectors_2010$ID
# for (detector in detectors_2010$ID) {
#   index <- which(det_bearings == detector, arr.ind = TRUE)
#   index[, 2] <- index[, 2] + 3
#   von_mises[index[, 1], detector] <- as.matrix(det_bearings)[index]
# }
# rm(index)
# 
# SE_bearings <- matrix(NA, nrow = nrow(det_bearings), ncol = nrow(detectors_2010))
# colnames(SE_bearings) <- detectors_2010$ID
# for (detector in detectors_2010$ID) {
#   index <- which(det_bearings == detector, arr.ind = TRUE)
#   index[, 2] <- index[, 2] + 2
#   SE_bearings[index[, 1], detector] <- as.matrix(det_bearings)[index]
# }
# rm(index)
# 
# time_of_arrivals <- matrix(NA, nrow = nrow(selection), ncol = nrow(detectors_2010))
# colnames(time_of_arrivals) <- detectors_2010$ID
# for (detector in detectors_2010$ID) {
#   index <- which(selection == detector, arr.ind = TRUE)
#   index[, 2] <- index[, 2] + 2
#   time_of_arrivals[index[, 1], detector] <- as.matrix(selection)[index]
# }
# rm(index)
# 
# sound_levels <- matrix(NA, nrow = nrow(selection), ncol = nrow(detectors_2010))
# colnames(sound_levels) <- detectors_2010$ID
# for (detector in detectors_2010$ID) {
#   index <- which(selection == detector, arr.ind = TRUE)
#   index[, 2] <- index[, 2] + 3
#   sound_levels[index[, 1], detector] <- as.matrix(selection)[index]
# }
# rm(index)
# 
# write_csv(as.data.frame(bearings_all), "Data/[PT] bearings_history_S5Y10_ALL.csv")
# # write_csv(as.data.frame(bearings_successful), "Data/bearings_history_s04y13_SUCCESSFUL.csv")
# # write_csv(as.data.frame(bearings_weight_1), "Data/bearings_history_s04y13_weight=1.csv")
# write_csv(as.data.frame(von_mises), "Data/[PT] von_Mises_history_S5Y10.csv")
# write_csv(as.data.frame(SE_bearings), "Data/[PT] se_bearings_history_S5Y10.csv")
# write_csv(as.data.frame(W), "Data/[PT] detection_history_S5Y10.csv")
# # write_csv(as.data.frame(W[dets$succesful_location == "successful", ]), 
# #           "Data/detection_history_s04y13_successful_location.csv")
# # write_csv(as.data.frame(W[dets$succesful_location != "successful", ]), 
# #           "Data/detection_history_s04y13_failed_location.csv")
# write_csv(as.data.frame(weights), "Data/[PT] weights_history_S5Y10.csv")
# write_csv(as.data.frame(sound_levels), "Data/[PT] sound_level_history_S5Y10.csv")
# write_csv(as.data.frame(time_of_arrivals), "Data/[PT] time_of_arrival_history_S5Y10.csv")
