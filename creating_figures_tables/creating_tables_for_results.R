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

df <- read.csv("Real data output/summary_estimates_bootstraps_model_33.csv")
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
df <- read.csv("Real data output/summary_real_bootstraps_model_33.csv")
output <- base::subset(df, select = c("point", "sd", "cv"))
row.names(output) <- df$parameter
table <- kable(x = output, 
               format = "latex",
               digits = 3,
               col.names = c("Point Estimate", "Standard Deviation", "Coefficient of Variation (%)"),
               label = "tab:best_model_results",
               caption = "Estimates for all quantities of interested for the best model, with Monte Carlo uncertainty estimates.")
print(table)

## -----------------------------------------------------------------------------
## Create a table of the results from the bootstraps for the supp. materials
## -----------------------------------------------------------------------------

load("Real data output/fits_1_35_nlminb_n=443.RData")
fit_overview <- data.frame(ID = 1:length(fits), 
                           density_model = sapply(fits, function(fit) Reduce(paste, deparse(fit$dat$f_density))),
                           N_hat = sapply(fits, function(fit) fit$N_total),
                           convergence = sapply(fits, function(fit) fit$fit_res$convergence),
                           n_par = sapply(fits, function(fit) length(fit$fit_res$par)),
                           llk = sapply(fits, function(fit) -fit$fit_res$objective),
                           AIC = sapply(fits, function(fit) fit$aic),
                           AICc = sapply(fits, function(fit) fit$aicc),
                           BIC = sapply(fits, function(fit) fit$bic))

fit_overview <- fit_overview[order(fit_overview$AIC), ]
output <- base::subset(fit_overview, select = c("ID", "density_model", "N_hat", 
                                                "n_par", "AIC"))
table <- kable(x = output, 
               format = "latex",
               digits = 3,
               col.names = c("ID", "Density Model", "N_hat", "#Par", "AIC"),
               label = "tab:results-from-the-999-bootstraps", row.names = FALSE)
## Print table to copy to LaTeX document on Overleaf
print(table)


### And now the parameter estimates
df <- read.csv("Real data output/summary_estimates_bootstraps_model_33.csv")
output_e <- base::subset(df, select = -c(8))
output_e$link <- c("iden.", "logit", rep("log", 6), "logit", rep("log", 11))
output_e <- output_e[, c(1, 9, 2:8)]

df <- read.csv("Real data output/summary_real_bootstraps_model_33.csv")
output_r <- base::subset(df, select = -c(8))
output_r$link <- c("iden.", "logit", rep("log", 6), "logit")
output_r <- output_r[, c(1, 9, 2:8)]

# output <- rbind(output_r, output_e[-c(1:9), ])
output <- output_e

output_sig <- output
output_sig[, c(3:9)] <- signif(output_sig[, c(3:9)], 2)

table <- kable(x = output, 
               format = "latex",
               digits = 5,
               col.names = c("Par.", "Link", "Est.", "LCL", "UCL", "LRCL", "URCL", 
                             "SD", "CV(%)"),
               label = "tab:best_model_overview_results")
## Print table to copy to LaTeX document on Overleaf
print(table)


df <- read.csv("Real data output/overview_estimates_bootstraps_model_33.csv")
