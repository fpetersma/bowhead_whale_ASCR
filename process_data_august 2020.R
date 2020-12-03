<<<<<<< Updated upstream
################################################################################
####  "process_data_august 2020.R"                                          ####
####  --------------------------------------------------------------------  ####
####  Turn new data from 27/08/2020 to a usable format, then writes it to   ####
####  csv files.                                                            ####
####  --------------------------------------------------------------------  ####
####  By: Felix Petersma                                                    ####
####  Last update: 27/08/2020                                               ####
################################################################################

#### Load libraries ############################################################

library(tidyverse)

#### Read in data ##############################################################
=======
# ============================================================================ #
# "process_data_august 2020.R"                                                 #
# ---------------------------------------------------------------------------- #
# Turn new data from 27/08/2020 to a usable format, then writes it to          #
# csv files.                                                                   #
# ---------------------------------------------------------------------------- #
# By: Felix Petersma                                                           #
# Last update: 27/08/2020                                                      #
# ============================================================================ #

# Load libraries ===============================================================

library(tidyverse)

# Define constants =============================================================

# choose whether calls should be filtered by successful localisation
USE_ONLY_LOC_CALLS <- FALSE 

# Read in data =================================================================
>>>>>>> Stashed changes

df <- read_csv("../JABES paper/Data/Data August 2020/Site5_Automated_2010.1_2minDASAR_250kmMaxRange_InfmUncertainty_20to300Hz_TL_PowerLaw15dB_20200826T151454.txt",
               col_names = c("n_sensors", 
                             "date_time_origin", 
                             "UTM_easting_origin",
                             "UTM_northing_origin",
                             "RMS_source_level", 
                             "RMS_discrepancy",
                             paste0(c("DASAR",
                                      "date_time_detection",
                                      "bearing",
                                      "range",
                                      "RMS_received",
                                      "RMS_noise",
                                      "peak_noise",
                                      "min_frequency_noise", 
                                      "max_frequency_noise"), 
                                    "_", 
                                    rep(1:6, each = 9)),
                             "NA"))

<<<<<<< Updated upstream
#### Start modifying data ######################################################

## Extract date and time from date-time
df_full <- df %>% 
  # mutate(date_time_origin = replace(date_time_origin, 
  #                                   date_time_origin == no_date_time, 
  #                                   NA)) %>%
  # select(1:date_time_origin, date_origin, everything()) %>% # reorder
  select(-61) %>% 
  mutate(successful = !is.nan(UTM_easting_origin),
         date = as.character(as.Date(date_time_origin))) %>% 
  select(1:date_time_origin, successful, everything()) # reorder

# unique(df_full$date)  

## If you want to only use successful location, run line below
df_full <- filter(df_full, successful)  %>% na_if(-1)

## Select 2010-08-31
df_full_31082010 <- filter(df_full, date == "2010-08-31")

## Extract the information
bearings <- select(df_full_31082010, paste0("bearing_", 1:6)) # %>% na_if(-1)

detections <- (!is.na(bearings)) * 1 # times 1 is a trick to make numeric

received_levels <- select(df_full_31082010, paste0("RMS_received_", 1:6)) %>% na_if(-1)

noise <- select(df_full_31082010, paste0("RMS_noise_", 1:6))

snr <- received_levels - noise

#### Write the data to .csv files ##############################################

write.csv(bearings, "../JABES paper/Data/bearings_31-08-2020_successful-loc.csv", row.names = F)
write.csv(detections, "../JABES paper/Data/detections_31-08-2020_successful-loc.csv", row.names = F)
write.csv(received_levels, "../JABES paper/Data/received_levels_31-08-2020_successful-loc.csv", row.names = F)
write.csv(noise, "../JABES paper/Data/noise_31-08-2020_successful-loc.csv", row.names = F)
write.csv(snr, "../JABES paper/Data/snr_31-08-2020_successful-loc.csv", row.names = F)

write.csv(df_full, "../JABES paper/Data/all-data_successful-loc.csv", row.names = F)
=======
# Start modifying data =========================================================

# extract date and time from date-time
df_full <- df %>% 
  select(-61) %>% 
  mutate(successful = !is.nan(UTM_easting_origin),
         date = as.character(as.Date(date_time_origin))) %>% 
  select(1:date_time_origin, successful, everything()) %>% #reorder
  slice(-c(1951,
           7118,
           15538,
           22025,
           40750,
           57241,
           69415,
           85813)) # remove faulty recordings

# unique(df_full$date)  

# if you want to only use calls with successful localisations, run line below
if (USE_ONLY_LOC_CALLS) {
  df_full <- filter(df_full, successful)  %>% na_if(-1)
}

# select the day of interest (2010-08-31)
df_full_31082010 <- filter(df_full, date == "2010-08-31")

# extract the relevant information
bearings <- select(df_full_31082010, paste0("bearing_", 1:6)) # %>% na_if(-1)
detections <- (!is.na(bearings)) * 1 # times 1 is a trick to make numeric
received_levels <- select(df_full_31082010, paste0("RMS_received_", 1:6)) %>% na_if(-1)
noise <- select(df_full_31082010, paste0("RMS_noise_", 1:6))
snr <- received_levels - noise

# Write the data to .csv files (remember to specify names correctly!) ==========

write.csv(bearings, "../JABES paper/Data/bearings_31-08-2020_all.csv", row.names = F)
write.csv(detections, "../JABES paper/Data/detections_31-08-2020_all.csv", row.names = F)
write.csv(received_levels, "../JABES paper/Data/received_levels_31-08-2020_all.csv", row.names = F)
write.csv(noise, "../JABES paper/Data/noise_31-08-2020_all.csv", row.names = F)
write.csv(df_full, "../JABES paper/Data/all-data_all.csv", row.names = F)
>>>>>>> Stashed changes
