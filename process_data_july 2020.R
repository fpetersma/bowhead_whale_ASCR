####### Turn new data from 24/07/2020 to a usable format ###########################

library(readr)
library(dplyr)
# df <- read.csv("Data/Data July 2020/Site5_Automated_2010.1_2minDASAR_250kmMaxRange_InfmUncertainty_20to300Hz_TL_PowerLaw15dB_20200723T174419.csv",
#                stringsAsFactors = F, header = F)
df <- read_csv("Data/Data July 2020/Site5_Automated_2010.1_2minDASAR_250kmMaxRange_InfmUncertainty_20to300Hz_TL_PowerLaw15dB_20200724T073309.csv",
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

unique(df_full$date)  
# df_full %>% select(date) %>% summarise(sum = sum(is.na(.)))

# # Get the dates for every call
# dates <- apply(df, 1, function(x) {
#   x <- x[paste0("date_time_detection_", 1:6)]
#   unique(as.character(as.Date(x[x != 0])))
# })
# df_full$date_origin <- head(dates, -1)

## If you want to only use successful location, run line below
df_full <- filter(df_full, successful)

## Select a time period
df_full <- filter(df_full, date == "2010-08-31")

## Extract the information
# Get the bearings
bearings <- select(df_full, paste0("bearing_", 1:6)) %>% na_if(-1)

detections <- (!is.na(bearings)) * 1 # times 1 is a trick to make numeric

received_levels <- select(df_full, paste0("RMS_received_", 1:6)) %>% na_if(-1)

noise <- select(df_full, paste0("RMS_noise_", 1:6))

snr <- received_levels - noise

plot(density(as.matrix(snr), na.rm = TRUE))

######## Exploring the range of source levels
density(df_full$RMS_source_level)


####### Turn data from 20/07/2020 to a usable format ###########################

library(readr)

df <- read_csv("Data/Data July 2020/old/Site5_Automated_2010.1_2minDASAR_250kmMaxRange_InfmUncertainty_20to300Hz_TL_PowerLaw15dB_20200719T090524.csv",
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

## Extract date and time from date-time

no_date_time <- df$date_time_origin[98232]

df_full <- df %>% 
  # mutate(date_time_origin = replace(date_time_origin, 
  #                                   date_time_origin == no_date_time, 
  #                                   NA)) %>%
  # select(1:date_time_origin, date_origin, everything()) %>% # reorder
  mutate(successful = date_time_origin != no_date_time) %>% 
  select(1:date_time_origin, successful, everything()) %>% # reorder
  head(-1) # remove faulty last row

# Get the dates for every call
dates <- apply(df, 1, function(x) {
  x <- x[paste0("date_time_detection_", 1:6)]
  unique(as.character(as.Date(x[x != 0])))
})
df_full$date_origin <- head(dates, -1)

## If you want to only use successful location, run line below
df_full <- filter(df_full, successful)

## Select a time period
df_full <- filter(df_full, date_origin == "2010-08-31")

detection_histories <- matrix(ncol = 6, 
                              nrow = nrow(df_full), 
                              dimnames = list(rownames(df_full), paste0("DASAR_", 1:6)))




