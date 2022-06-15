## =============================================================================
## Creating plots for parameter estimates, and other results presented in the 
## manuscript. 
##
## First created on 21/02/2021 by Felix T Petersma
## =============================================================================
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(ggtext)

## -----------------------------------------------------------------------------
## A concept of what a plot of the parameter estimates for the case study could
## look like. A lot of it is based on
##    https://ggplot2.tidyverse.org/reference/geom_linerange.html
## with slight modification.
## 
## I am satisfied with the looks of it, somewhat at least.
## [23/02/22] Len suggested doing violin plots, since those convey more info.
##            I added a violin plot version.
## -----------------------------------------------------------------------------

df <- read.csv("Real data output/summary_estimates_bootstraps_model_31.csv")

## Create a crossbar
ggplot(data = df, 
       aes(y = factor(parameter, levels = rev(parameter)), # specify levels this way to retain the correct order
           x = rep(0, 17))) + # plot 17 zeros as place holders
  geom_crossbar(aes(xmax = rel_97.5. * 100, xmin = rel_2.5. * 100), 
                colour = "azure4", fill = "azure3", alpha = 0.7, size = 0.5,
                fatten = NULL) + ## Use fatten = NULL for no mean/middle line
  theme_light() +
  geom_vline(xintercept = 0, colour = "red", linetype = "dashed") + # not sure if I want this line
  theme(legend.position = "none") +
  labs(y = "Estimated parameter", x = "95% confidence interval, relative to point estimate")

## Create a violin plot version
df <- read.csv("Real data output/overview_rel_estimates_bootstraps_model_31.csv")
# Remove 


df_melt <- melt(df)

ggplot(data = df_melt) +
  geom_violin(aes(x = variable, y = value)) +
  theme_light() + 
  coord_flip()
  geom_vline(xintercept = 0, colour = "red", linetype = "dashed") + # not sure if I want this line
  theme(legend.position = "none") +
  labs(y = "Estimated parameter", x = "95% confidence interval, relative to point estimate")


## -----------------------------------------------------------------------------
## The plot below is an idea for plotting the relative error (RE) and coefficient
## of variation (CV).
##
## The RB will be a range, as it can be derived for all 100 simulations. 
## 
## The CV will be a point estimate, as it requires the Monte Carlo standard
## deviation. 
## -----------------------------------------------------------------------------

all_N <- matrix(NA, nrow = 100, ncol = 12)
all_N <- as.data.frame(all_N)

## VARIABLE SOURCE LEVEL SIMULATIONS -------------------------------------------
## Add the true realised total abundance (N)
load("1000_simulations_fixedSL=FALSE.RData")
all_N[, 1] <- sapply(simulated_data[1:100], function(dat) dat$N_expected)
all_N[, 2] <- sapply(simulated_data[1:100], function(dat) dat$N_realised)

## Variable source level fits
load("Simulation output/fits_with_fixedSL=FALSE_USE_BEARINGS=2_1-100_on_simulations_with_fixedSL=FALSE.RData")
all_N[, 3] <- sapply(results, function(fit) fit$N_total)

## Fixed source level fits
load("Simulation output/fits_with_fixedSL=TRUE_USE_BEARINGS=2_1-100_on_simulations_with_fixedSL=FALSE.RData")
all_N[, 4] <- sapply(results, function(fit) fit$N_total)

## Simple bearings model
load("Simulation output/fits_with_fixedSL=FALSE_USE_BEARINGS=1_1-100_on_simulations_with_fixedSL=FALSE.RData")
all_N[, 5] <- sapply(results, function(fit) fit$N_total)

## No bearings included
load("Simulation output/fits_with_fixedSL=FALSE_USE_BEARINGS=0_1-100_on_simulations_with_fixedSL=FALSE.RData")
all_N[, 6] <- sapply(results, function(fit) fit$N_total)

## Constant density model
load("Simulation output/fits_with_fixedSL=FALSE_USE_BEARINGS=2_constant_density_1-100_on_simulations_with_fixedSL=FALSE.RData")
all_N[, 7] <- sapply(results, function(fit) fit$N_total)

## FIXED SOURCE LEVEL SIMULATIONS ----------------------------------------------
## Add the true realised total abundance (N)
load("1000_simulations_fixedSL=TRUE.RData")
all_N[, 8] <- sapply(simulated_data[1:100], function(dat) dat$N_expected)
all_N[, 9] <- sapply(simulated_data[1:100], function(dat) dat$N_realised)

## Fixed source level fits
load("Simulation output/fits_with_fixedSL=TRUE_USE_BEARINGS=2_1-100_on_simulations_with_fixedSL=TRUE.RData")
all_N[, 10] <- sapply(results, function(fit) fit$N_total)

## Variable source level fits
load("Simulation output/fits_with_fixedSL=FALSE_USE_BEARINGS=2_1-100_on_simulations_with_fixedSL=TRUE_density_pars=c(-16,57,-68.5).RData")
all_N[, 11] <- sapply(results, function(fit) fit$N_total)

## Simple bearings model
load("Simulation output/fits_with_fixedSL=TRUE_USE_BEARINGS=1_1-100_on_simulations_with_fixedSL=TRUE.RData")
all_N[, 12] <- sapply(results, function(fit) fit$N_total)

## No bearings included
load("Simulation output/fits_with_fixedSL=TRUE_USE_BEARINGS=0_1-100_on_simulations_with_fixedSL=TRUE.RData")
all_N[, 13] <- sapply(results, function(fit) fit$N_total)

## Constant density model
load("Simulation output/fits_with_fixedSL=TRUE_USE_BEARINGS=2_constant_density_1-100_on_simulations_with_fixedSL=TRUE.RData")
all_N[, 14] <- sapply(results, function(fit) fit$N_total)

## Make matrix pretty and clear environment
colnames(all_N) <- c("N_expected_var_SL",
                     "N_realised_var_SL",
                     "scenario_1.1",
                     "scenario_1.2",
                     "scenario_1.3",
                     "scenario_1.4",
                     "scenario_1.5",
                     "N_expected_fixed_SL",
                     "N_realised_fixed_SL",
                     "scenario_2.1",
                     "scenario_2.2",
                     "scenario_2.3",
                     "scenario_2.4",
                     "scenario_2.5")

rm(simulated_data, results)

## EXTRACT RELATIVE BIAS AND COEFFIENT OF VARIATION ----------------------------
## Do we use the N_expected or N_realised as the 'truth' in the bias derivation?
## Current version uses the expectation.
RE <- all_N
RE[, 3:7] <- (RE[, 3:7] - RE[, 1]) / RE[, 1] * 100
RE[, 10:14] <- (RE[, 10:14] - RE[, 8]) / RE[, 8] * 100
RE <- RE[, -c(1, 2, 8, 9)]

CV <- data.frame(CV = apply(all_N[, -c(1, 2, 8, 9)], 2, function(x) {
  return(sd(x) / mean(x))
}))

RB <- data.frame(relative_bias = colMeans(RE),
                 t(apply(RE, 2, function(x) {
                   quantile(x, probs = c(0.025, 0.975))
                 })))
colnames(RB) <- c("relative_bias", "lower_95", "upper_95")

## Write the files to csv files, so that other scripts can use them
write.csv(x = all_N, file = "Simulation output/overview_N_simulation_scenarios.csv",
          row.names = FALSE)
write.csv(x = RE, file = "Simulation output/summary_relative_error.csv",
          row.names = FALSE)
write.csv(x = t(CV), file = "Simulation output/summary_coefficient_of_variation.csv",
          row.names = FALSE)
write.csv(x = RB, file = "Simulation output/summary_relative_bias.csv",
          row.names = FALSE)

## CREATE THE PLOTS ------------------------------------------------------------
## Prepare the data
df <- cbind(t(RE), CV, RB = RB$relative_bias, scenario = rownames(t(RE)))
df_melt <- melt(df, 
                value.name = "RE", 
                #measure.vars = row.names(RE),
                id = c("CV", "RB", "scenario"))

## Create black and white boxplots
p_bw <- ggplot(data = df_melt, mapping =  aes(x = scenario, y = RE)) +
  geom_violin(trim = FALSE, width = 1.5, lwd = 1,
              colour = "grey75", 
               fill = "grey75", 
               alpha = 1) +
  stat_summary(fun = mean, geom = "point", shape = 16) +
  labs(y = expression(paste("Relative error in ", hat(N)," (%)")), 
       x = "Simulation scenario and model used for analysis") +
  geom_hline(yintercept = 0, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 5.5, colour = "grey", linetype = "79") +
  # theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
  scale_x_discrete(labels = c("1a", "1b", "1c", "1d", "1e", 
                              "2a", "2b", "2c", "2d", "2e")) +
  scale_y_continuous(breaks = seq(from = -300, to = 900, by = 100),
                     labels = c("CV", "RB", seq(from = -100, to = 900, by = 100))) +
  theme_light() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks = element_blank()) + #,
  # axis.ticks.y = element_line(colour = c("transparent", rep("black", 11)))) +
  # geom_label(mapping = aes(y = -200, label = paste0(round(CV * 100, 0), "%")),
  #                         fill = "white",
  #                         size = 3.5,
  #           alpha = 0.001)
  geom_hline(yintercept = c(-200, -300), colour = "white") +
  annotate(geom = "label", x = 1:10, y = -300,
           label = paste0(round(CV$CV * 100, 0), "%"),
           fill = "white", label.size = NA) +
  annotate(geom = "label", x = 1:10, y = -200,
           label = paste0(round(RB$relative_bias, 0), "%"),
           fill = "white", label.size = NA) +
  annotate(geom= "label", x = c(3, 8), y = 900, 
           label = c("Variable Source Level", "Fixed Source Level"), 
           fill = "grey", label.size = NA, colour = "white")
p_bw
# +
  # geom_text(aes(x = variable, y = -150, label = )) 
  # geom_point(aes(x = variable, y = -150))


## Create the boxplots plots
p_colour <- ggplot(data = RE_melt) +
  geom_boxplot(mapping = aes(x = variable, y = RE),
               colour = "azure4", 
               fill = "azure3", 
               alpha = 0.8,
               size = 0.5) +
  labs(y = "Relative error (%)", 
       x = "Simulation scenario") +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
  # theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
  scale_x_discrete(labels = c("1.1", "1.2", "1.3", "1.4", "1.5", 
                              "2.1", "2.2", "2.3", "2.4", "2.5")) +
  scale_y_continuous(breaks = seq(from = -100, to = 900, by = 100)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank()) +
  geom_point(aes(x = variable, y = -150))

p2 <- ggplot() +
  geom_point(aes(x = unique(RE_melt$variable), y = -150))
p2


## =============================================================================
## Making a plot of the detection function and the source level distribution.
## =============================================================================

## Source level
## -----------------------------------------------------------------------------

## Set parameters
mu_s <- 158.838
sd_s <- 5.419

x_sl <- 135:185
y_sl <- dnorm(x_sl, mu_s, sd_s)

p_sl <- ggplot(mapping = aes(x = x_sl, y = y_sl)) +
  geom_area(fill = "grey", alpha = 0.6) +
  geom_line(colour = "grey", lwd = 1) +
  theme_classic() + #theme_hc() +
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank()
        ) +
  labs(x = "Source level (dB)", y = "Density") +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(expand = c(0, 0))

## Detection function 
## -----------------------------------------------------------------------------

g0 <- 0.546
trunc_level <- 96

x <- c(80:95, 95.99, 96:110)
y <- ifelse(x < trunc_level, 0, g0)

p_df <- ggplot(mapping = aes(x = x, y = y)) +
  geom_step(colour = "grey", lwd = 1.5) +
  theme_classic() + 
  # theme(axis.ticks.y = element_blank()) +
  labs(x = "Received level (dB)", y = "Detection probability") +
  geom_area(fill = "grey", alpha = 0.6) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
  scale_x_continuous(expand = c(0, 0), limits = c(90, 105))

ggarrange(p_sl, p_df, ncol = 2, nrow = 1, labels = c("A", "B"))


