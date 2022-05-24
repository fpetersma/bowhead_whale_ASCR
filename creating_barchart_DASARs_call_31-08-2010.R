## Creating a barchart of DASARs invovled per call for 31-08-2010

## LOAD SINGLETON DATA


library(readr)
library(dplyr)

########################## Load the .tsv files #################################
folder_loc <- "../Export all Aaron/Export_Automated_2010_Initial/2010 [no site 2]/" 
files <- list.files(folder_loc)

site <- 5 # Change this for different sites
n_detectors <- 6

index <- substr(files, 3, 8) == "0831S5"

files <- files[index] # Only keep files belonging to specific 

# only keep singletons
files <- files[2]

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

singletons <- detections

singletons_above_tr <- singletons[singletons$rms_db_re_1mPa_1 >= 96, ]

single_dets <- rep(1, nrow(singletons_above_tr))

## So far, we have the original singletons. However, we also have the newly created
## singletons once we truncate the data on 96dB received level. 

received_levels <- read.csv("Data/received_levels_31-08-2010_all.csv")

sufficient_rl <- received_levels  >= 96 # get indices 
sufficient_rl[is.na(sufficient_rl)] <- FALSE

dets <- sufficient_rl * 1
enough_dets <- Rfast::rowsums(dets) > 0
dets <- dets[enough_dets, ]

multi_dets <- Rfast::rowsums(dets)

table(single_dets)
table(multi_dets)

dets_involved <- as.data.frame(table(c(single_dets, multi_dets)))

##### THIS IS THE POINT WHERE I WANTED TO GET:
### 1 - 1493
### 2 - 284
### 3 - 111
### 4 - 37
### 5 - 23
### 6 - 5

dets_involved$percentages <- round(dets_involved$Freq / 
                                     sum(dets_involved$Freq), 3) * 100
dets_involved$n_DASAR <- 1:6

# hist(rowsums(as.matrix(W)))
ggplot(data = dets_involved, mapping =  aes(x = n_DASAR, y = Freq)) + 
  theme_classic() +
  geom_bar(stat = "identity", fill = "grey75") +
  labs(x = "DASARs involved", y = "Number of calls") +
  scale_x_continuous(breaks = 1:6) +
  geom_text(aes(label = paste0(percentages, "%")), 
            position=position_dodge(width = 0.9), vjust = -0.25) +
  geom_line(data = data.frame(x = c(0.55, 1.45), y = c(587, 587)),
            mapping = aes(x = x, y = y), linetype = "dashed", colour = "black") +
  ylim(c(0,1600))



