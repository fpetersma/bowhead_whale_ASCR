################################################################################
####  "source_level_distribution.R"                                         ####
####  --------------------------------------------------------------------  ####
####  Explore the distribution of source levels. These source levels were   ####
####  derived by Aaron Thode using an estimate of the location and an       ####
####  estimate of the propagation loss.                                     ####
####  --------------------------------------------------------------------  ####
####  By: Felix Petersma                                                    ####
####  When: 25/08/2020                                                      ####
################################################################################

#### Load libraries ############################################################

library(tidyverse)
library(ggsci)
library(MASS)

#### Read in data ##############################################################

df_all <- read_csv("../JABES paper/Data/all-data_successful-loc.csv")

#### Analyse source level distribution #########################################

# Plot of all data by number of detections
ggplot(data = df_all) + 
  geom_histogram(mapping = aes(x = RMS_source_level, group = n_sensors,
                             fill = as.factor(n_sensors)), alpha = 1, bins = 200) + 
  theme_bw() +
  ggtitle("Histograms of source levels for site 5, 2010.", 
          subtitle = "Grouped by number of involved DASARs.") +
  scale_fill_npg(name = "Number of DASARs \ninvolved") + labs(y = "Number of calls", x = "Source level (RMS) in db re 1microPa") 
summary(df_all$RMS_source_level) 

quantile(df_all$RMS_source_level, probs = c(0.005, 0.995)) 

fitdistr(df_all$RMS_source_level, densfun = "normal")

# Plot of the data for 2010-08-31
df_all %>% 
  filter(date == "2010-08-31") %>% 
  ggplot() +
  geom_density(mapping = aes(x = RMS_source_level, group = n_sensors,
                             fill = n_sensors), alpha = 0.5) + theme_minimal()
df_all %>%
  filter(date == "2010-08-31") %>%
  select(RMS_source_level) %>% 
  summary()

quantile(df_all$RMS_source_level[df_all$date == "2010-08-31"],
         probs = c(0.005, 0.995))

fitdistr(df_all$RMS_source_level[df_all$date == "2010-08-31"], densfun = "normal")

#### Conclusion ################################################################
# In conclusion, it seems that the source levels vary vastly, with a very 
# similar distribution overall as just for 31-08-2010. In both cases, the 
# smallest intervals that contain at least 99% of the data are roughly c(126, 
# 175). This means that, if I would use c(120, 180) as bounds for my 
# integration, I will have at least 99% of the data. We use a similar method
# to check whether the integration area in space is large enough; there, we 
# check whether the probability of detection is smaller than 1% at the bounds.
# The fitted normal distributions are N(153.7, 8.1) and N(153.1, 7.7), 
# respectively.
  
