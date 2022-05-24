## ========================================================================== ##
## creating_tables_for_results.R                                              ##
##                                                                            ##
## In this script I will build tables that can be imported in LaTeX.          ##
## ========================================================================== ##

library(tidyverse)
library(knitr)

## -----------------------------------------------------------------------------
## Create a table with relative bias (= mean relative error) and empirical cv
## for all 10 simulation scenarios. 
## -----------------------------------------------------------------------------

## Load the data
df <- read.csv("Simulation output/overview_N_simulation_scenarios.csv")
df_means <- colMeans(df)
RB <- df_means

## Derive the relative bias 
RB[3:7] <- (RB[3:7] - RB[1]) / RB[1] * 100
RB[10:14] <- (RB[10:14] - RB[8]) / RB[8] * 100
RB <- as.data.frame(t(RB[-c(1, 2, 8, 9)]))

## Load CV estimates
CV <- read.csv("Simulation output/summary_coefficient_of_variation.csv")

## Create table
joint <- t(rbind(RB, CV))
colnames(joint) <- c("RB", "CV")
table <- kable(x = joint, 
               format = "latex",
               digits = 3,
               col.names = c("Relative Bias (%)", "Coefficient of Variation (%)"),
               label = "tab:simulation_rb_cv",
               caption = "Uncertainty around hat{N} from 100 fits to simulated data in 10 scenarios.")

## Print table to copy to LaTeX document on Overleaf
print(table)

## -----------------------------------------------------------------------------
## Create a table of the results from the best model, and the Monte Carlo
## uncertainty.
## -----------------------------------------------------------------------------

df <- read.csv("Real data output/summary_estimates_bootstraps_model_31.csv")
output <- base::subset(df, select = c("point", "sd", "cv"))
row.names(output) <- df$parameter
table <- kable(x = output, 
               format = "latex",
               digits = 3,
               col.names = c("Point Estimate", "Standard Deviation", "Coefficient of Variation (%)"),
               label = "tab:best_model_results",
               caption = "Estimates for all quantities of interested for the best model, with Monte Carlo uncertainty estimates.")

## Print table to copy to LaTeX document on Overleaf
print(table)


## Do the same for the real data
df <- read.csv("Real data output/summary_real_bootstraps_model_31.csv")
output <- base::subset(df, select = c("point", "sd", "cv"))
row.names(output) <- df$parameter
table <- kable(x = output, 
               format = "latex",
               digits = 3,
               col.names = c("Point Estimate", "Standard Deviation", "Coefficient of Variation (%)"),
               label = "tab:best_model_results",
               caption = "Estimates for all quantities of interested for the best model, with Monte Carlo uncertainty estimates.")
print(table)
