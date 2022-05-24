## ========================================================================== ##
## In this script, I extract all relevant info from bootstraps on model 31.   ##
## These extracts are then stored in .csv files which are easier to work with.##                                                                 ##
## ========================================================================== ##

## Turn off scientific notation
options(scipen=999)

## -----------------------------------------------------------------------------
## 1. Load the data and prepare for extraction
## -----------------------------------------------------------------------------

## Load model 31
load("Real data output/33_fits_real_data_SL=100-220.RData")
best_model <- fits[[31]]

rm(fits)

## Load the 999 fits on bootstrapped data
load("Real data output/bootstraps on best model (model 31) - unfinished/results_1-999.RData")

total <- total_results # rename
rm(total_results)

## Add the best model in position 1000
total[[1000]] <- best_model

## The bootstraps are missing the mesh, so add that for completeness
for (i in seq(999)) {
  x <- total[[i]]
  x$mesh <- total[[1000]]$mesh
  total[[i]] <- x
}

## Clean the environment
rm(i, x, best_model)

## -----------------------------------------------------------------------------
## 2. Extract parameter and total abundance estimates with confidence intervals
## -----------------------------------------------------------------------------

## Create data frame to store all parameter estimates in
df_estimates <- data.frame(N = NA,
                           logit_g0 = NA,
                           log_beta_r = NA,
                           log_sd_r = NA,
                           log_mu_s = NA,
                           log_sd_s = NA,
                           log_kappa_low = NA,
                           log_kappa_high = NA,
                           logit_mix_bear = NA,
                           beta_intercept = NA,
                           beta_depth = NA,
                           beta_depth2 = NA,
                           beta_sd2c_1 = NA,
                           beta_sd2c_2 = NA,
                           beta_sd2c_3 = NA,
                           beta_sd2c_4 = NA,
                           beta_sd2c_5 = NA)

## Loop through all bootstraps and the best model 
for (i in seq(1000)) {
  df_estimates[i, "N"] <- total[[i]]$N_total
  df_estimates[i, 2:ncol(df_estimates)] <- total[[i]]$fit_res$par
}

## Relative estimates, relative to point estimate of the best model
df_rel_estimates <- as.data.frame(apply(df_estimates, 2, function(x) {
  return((x - x[1000]) / abs(x[1000]))
}))

# ## Reverse row order to have the best model as the first observation
# df_estimates <- df_estimates[rev(row.names(df_estimates)), ]

## Create a summary data frame
summary_estimates <- data.frame(parameter = colnames(df_estimates), 
                                point = unlist(df_estimates[1000, ]))

## Add 95% Monte Carlo intervals
for (i in seq(nrow(summary_estimates))) {
  summary_estimates[i, 3:4] <- quantile(df_estimates[, i], c(0.025, 0.975))
}
colnames(summary_estimates)[3:4] <- names(quantile(df_estimates[, 1], 
                                                   c(0.025, 0.975)))

## Add relative 95% Monte Carlo intervals
for (i in seq(nrow(summary_estimates))) {
  summary_estimates[i, 5:6] <- quantile(df_rel_estimates[, i], c(0.025, 0.975))
}
colnames(summary_estimates)[5:6] <- paste0("rel_", names(quantile(df_estimates[, 1], 
                                                                  c(0.025, 0.975))))

## Add standard deviation, variance and cv
summary_estimates[["sd"]] <- apply(df_estimates, 2, sd)
summary_estimates[["variance"]] <- summary_estimates$sd^2
summary_estimates[["cv"]] <- summary_estimates$sd / abs(apply(df_estimates, 2, mean)) * 100

## -----------------------------------------------------------------------------
## 3. Now do the same for the real parameters
## -----------------------------------------------------------------------------
df_real <- data.frame(N = NA,
                      g0 = NA,
                      beta_r = NA,
                      sd_r = NA,
                      mu_s = NA,
                      sd_s = NA,
                      kappa_low = NA,
                      kappa_high = NA,
                      mix_bear = NA)

## Loop through all bootstraps and the best model 
for (i in seq(1000)) {
  df_real[i, "N"] <- total[[i]]$N_total
  df_real[i, 2:ncol(df_real)] <- total[[i]]$real_par[1:8]
}

## Relative real, relative to point estimate of the best model
df_rel_real <- as.data.frame(apply(df_real, 2, function(x) {
  return((x - x[1000]) / abs(x[1000]))
}))

# ## Reverse row order to have the best model as the first observation
# df_real <- df_real[rev(row.names(df_real)), ]

## Create a summary data frame
summary_real <- data.frame(parameter = colnames(df_real), 
                           point = unlist(df_real[1000, ]))

## Add 95% Monte Carlo intervals
for (i in seq(nrow(summary_real))) {
  summary_real[i, 3:4] <- quantile(df_real[, i], c(0.025, 0.975))
}
colnames(summary_real)[3:4] <- names(quantile(df_real[, 1], 
                                              c(0.025, 0.975)))

## Add relative 95% Monte Carlo intervals
for (i in seq(nrow(summary_real))) {
  summary_real[i, 5:6] <- quantile(df_rel_real[, i], c(0.025, 0.975))
}
colnames(summary_real)[5:6] <- paste0("rel_", names(quantile(df_real[, 1], 
                                                             c(0.025, 0.975))))

## Add standard deviation, variance and cv
summary_real[["sd"]] <- apply(df_real, 2, sd)
summary_real[["variance"]] <- summary_real$sd^2
summary_real[["cv"]] <- summary_real$sd / abs(apply(df_real, 2, mean)) * 100

## -----------------------------------------------------------------------------
## 4. Create csv files
## -----------------------------------------------------------------------------

write.csv(df_estimates, 
          file = "Real data output/overview_estimates_bootstraps_model_31.csv", 
          row.names = FALSE)

write.csv(df_rel_estimates, 
          file = "Real data output/overview_rel_estimates_bootstraps_model_31.csv", 
          row.names = FALSE)

write.csv(summary_estimates, 
          file = "Real data output/summary_estimates_bootstraps_model_31.csv", 
          row.names = FALSE)

write.csv(df_real, 
          file = "Real data output/overview_real_bootstraps_model_31.csv", 
          row.names = FALSE)

write.csv(df_rel_real, 
          file = "Real data output/overview_rel_real_bootstraps_model_31.csv", 
          row.names = FALSE)

write.csv(summary_real, 
          file = "Real data output/summary_real_bootstraps_model_31.csv", 
          row.names = FALSE)
