## =============================================================================
## WHAT DOES THIS SCRIPT DO?
## =============================================================================

## run the sourceCpp() for test runs
# Rcpp::sourceCpp("Scripts/bowhead whales/nllRcpp.cpp")

## =============================================================================
## 1. LOAD DATA AND LIBRARIES
## -----------------------------------------------------------------------------

## Libraries
library(Rcpp)       # for combining R and cpp
library(parallel)   # for parallel processing
library(mgcv)


## Load the real data
## DO LATER!


## Load RData file with data from main.R
load("../../TMB ft. Fanny-E/TMB_testing/test_data.RData")

## Scale the density predictors to improve convergence (think about this --
## does not apply to splines)
# gam_fit <- dat$gam_fit
design_matrix <- dat$design
design_matrix[, -1] <- scale(design_matrix[, -1], FALSE, TRUE)
# design_pars <- rep(0, ncol(design_matrix))
# names(design_pars) <- colnames(design_matrix)

mesh <- read.csv("Data/alaska_albers_grid_adaptive_levels=2_inner=10k_outer=50k_maxD2C=Inf_area=8450_n=438.csv")

## Create standardised covariates data.frame
covariates <- dplyr::select(mesh, c(distance_to_coast, depth)) %>% 
  mutate(distance_to_coast2 = distance_to_coast ^ 2, 
         depth = abs(depth),
         depth2 = depth ^ 2,
         logdepth = log(depth)) 
dat$covariates <- scale(covariates)

dat$A_x <- dat$A_x$area

## Constants
A_s <- 3
FIXED_SL <- TRUE
USE_BEARINGS <- 2

## Run for test data
# load("C:/Users/felix/OneDrive - University of St Andrews/Documents/University of St Andrews/PhD/TMB ft. Fanny-E/TMB_testing/simulated_test_data_tl=75_ngrid=440_density=-2_g0=0.5_mu_s=140_betar=15_sdr=3_kappa_low=1_kappa_high=19_mix_par=0.1_N=1091.RData")
# dat$det_hist <- dat_sim$det_hist
# dat$received_levels <- dat_sim$received_levels
# bearings_deg <- circular::circular(dat_sim$bearings,
#                           units = "degrees",
#                           template = "geographics") # 'geographics' means clockwise rotation with 0 at north
# dat$bearings_rad <- circular::conversion.circular(bearings_deg)
# dat$trunc_level <- 75

## =============================================================================
## RUN WHEN USING SIM_DATA
## -----------------------------------------------------------------------------
## Libraries
library(Rcpp)       # for combining R and cpp
library(parallel)   # for parallel processing
library(mgcv)

dat <- sim_data

f_density <- dat$f_density
## Constants
A_s <- 3
FIXED_SL <- sim_data$FIXED_SL
USE_BEARINGS <- 2

## =============================================================================


## =============================================================================
## SET UP DENSITY MODEL
## -----------------------------------------------------------------------------
# f_density <- D ~ distance_to_coast + distance_to_coast2
# f_density <- D ~ 1
# f_density <- dat$f_density
## Create GAM object with scaled covariates
gam_fit <- gam(formula = f_density,
               data = data.frame(D = 0, dat$covariates))

## Create data object
data_obj <- list(Y_rec = dat$bearings_rad,      # in radians
                 Y_grid = dat$grid_bearings,    # in radians
                 X = dat$distances,
                 W = dat$det_hist,
                 R = dat$received_levels,
                 A_x = dat$A_x,
                 f_density = f_density,
                 # A_s = A_s,
                 trunc_level = dat$trunc_level,
                 # design_matrix = design_matrix,
                 gam_fit = gam_fit,
                 # S = seq(from = 80 + A_s / 2, to = 200 - A_s / 2, by = A_s),
                 n_det = ncol(dat$distances),
                 n_grid = nrow(dat$distances),
                 n_call = nrow(dat$det_hist),
                 FIXED_SL = FIXED_SL,
                 USE_BEARINGS = USE_BEARINGS) # create data
if (!FIXED_SL) {
  data_obj$A_s <- A_s
  # data_obj$S <- seq(from = 80 + A_s / 2, to = 200 - A_s / 2, by = A_s) # for midpoint rule
  data_obj$S <- seq(from = 80, to = 200, by = A_s)                     # for Simpson's rule
  data_obj$n_sl <- length(data_obj$S)
}
data_obj$min_det_matrix <- diag(1, data_obj$n_det + 1, data_obj$n_det) # matrix used for .p 

## Create start parameter value
density_pars <- rep(-1, ncol(model.matrix(data_obj$gam_fit)))
names(density_pars) <- colnames(model.matrix(data_obj$gam_fit))

parameters <- c(logit_g0 = 0 , 
                log_beta_r = log(15) , 
                log_sd_r = log(2) , 
                log_mu_s = log(150))
if (!FIXED_SL) parameters <- c(parameters,
                               log_sd_s = 2)
if (USE_BEARINGS == 1) parameters <- c(parameters,
                                       log_kappa = 1)
if (USE_BEARINGS == 2) parameters <- c(parameters,
                                       log_kappa_low = 1 ,
                                       log_kappa_high = 2 ,
                                       logit_mix_bear = -0.5)
parameters <- c(parameters, density_pars)

data_obj$dens_par_names <- colnames(model.matrix(data_obj$gam_fit))
## =============================================================================

## =============================================================================
## 2. SOURCE OTHER R AND CPP SCRIPTS
## -----------------------------------------------------------------------------
source("Scripts/bowhead whales/nllR.R")

## =============================================================================

## =============================================================================
## BONUS: BENCHMARK A SINGLE RUN
## -----------------------------------------------------------------------------
microbenchmark::microbenchmark({
  test <- nllR(parameters, data_obj)
}, times = 5)

## NOTES: 1 core(s) = 13.16 seconds; 
##        2 core(s) =  7.26 seconds; 
##        4 core(s) =  5.03 seconds. 

## =============================================================================
## 3. RUN THE OPTIMISER 
## -----------------------------------------------------------------------------
# res <- optim(par = parameters, fn = nllR, dat = data_obj, method = "BFGS",
#              control = list(trace = 6, REPORT = TRUE))
res <- nlminb(start = parameters, objective = nllR, dat = data_obj,
              lower = -50,
              upper = 50,
              control = list(trace = 1))

## =============================================================================
#  Truth:             0.405     2.701   0.693     5.011   0.000       2.944   -2.197   -1        2       -1.5
#  35:     4263.9473: 0.452523  2.71824 0.700802  5.01431 0.00945608  2.82123 -2.38982 -1.13384  2.48271 -2.22861
#  33:     4300.1762: 0.453393  2.72198 0.718492  5.01494 -0.198829   2.93911 -2.20946 -1.20934  2.86865 -2.69616
#  34:     4226.5505: 0.477780  2.72429 0.673404  5.01761 -0.177191   2.90168 -2.35719 -1.07825  2.10859 -1.71153
## =============================================================================
## 4. ADD OTHER INFORMATION TO DATA OBJECT
## -----------------------------------------------------------------------------
fit_obj <- list(dat = data_obj,
                start_pars = parameters,
                fit_res = res)
real_par <- exp(res$par)
real_par[which(names(res$par) %in% c("logit_g0", "logit_mix_par"))] <- 
  exp(res$par[which(names(res$par) %in% c("logit_g0", "logit_mix_par"))]) /
  (1 + exp(res$par[which(names(res$par) %in% c("logit_g0", "logit_mix_par"))]))
fit_obj$real_par <- real_par
fit_obj$aic <- -2 * res$value + 2 * length(res$par) # Burnham & Anderson (2002)
fit_obj$aicc <- fit_obj$aic + 2 * length(res$par) * (length(res$par) + 1) / 
  (data_obj$n_call - length(res$par) - 1) # Burnham & Anderson (2002)
fit_obj$bic <- -2 * res$value + length(res$par) * log(data_obj$n_call)

gam_fit <- data_obj$gam_fit
gam_fit$coefficients <- res$par[names(gam_fit$coefficients)]
fit_obj$D <- data.frame(area = data_obj$A_x, 
                        density = exp(predict(gam_fit, data = model.matrix(gam_fit))))
fit_obj$D$N <- fit_obj$D$area * fit_obj$D$density
fit_obj$N_total <- sum(fit_obj$D$N)

## =============================================================================