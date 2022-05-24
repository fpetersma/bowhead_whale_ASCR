library(Rcpp)

## Load RData file with data from main.R, and extract data and parameters for TMB
load("../../TMB ft. Fanny-E/TMB_testing/test_data.RData")

## Scale the density predictors to improve convergence (think about this -- does not apply to splines)
design_matrix <- dat$design
design_matrix[, -1] <- scale(design_matrix[, -1], TRUE, TRUE)
design_pars <- rep(0, ncol(design_matrix))
names(design_pars) <- colnames(design_matrix)

A_s <- 5

data <- list(Y_rec = dat$bearings_rad,
             Y_grid = dat$grid_bearings,
             X = dat$distances,
             W = dat$det_hist,
             R = dat$received_levels,
             A = dat$A_x$area,
             trunc_level = dat$trunc_level,
             design_matrix = dat$design,
             S = seq(from = 80 + A_s / 2, to = 200 - A_s / 2, by = A_s)) # create data

## Create start parameter value

density_pars <- rep(0, ncol(design_matrix))
names(density_pars) <- colnames(design_matrix)

parameters <- c(logit_g0 = 0.6 , 
                # log_kappa = 3, 
                log_beta_r = log(18) , 
                log_sd_r = 1 , 
                log_mu_s = log(160) ,
                log_sd_s = 2 ,
                log_kappa_low = 1 ,
                log_kappa_high = 2.0 ,
                logit_mix_bear = -0.4 ,
                density_pars) # define starting values for parameters


## Create CPP function
sourceCpp("Scripts/bowhead whales/llkRcpp.cpp")

## Create R wrapper
llkR <- function(par, dat) {
  print(par)
  pars <- as.list(par[1:8])
  pars[["density_pars"]] <- par[9:length(parameters)]
  
  nll <- llkRcpp(data = dat, par = pars)
  return(nll)
}
microbenchmark::microbenchmark({
  llkR(parameters, data)
}, times = 3)

rm("dat", "par")

microbenchmark::microbenchmark({
  opt <- optim(par = parameters, fn = llkR, method = "BFGS", 
               control = list(trace = 1, REPORT = 1), dat = data)
})



## Do TMB things
obj <- MakeADFun(data, parameters, DLL="ascr_old_marginalise", silent = FALSE)
# obj <- normalize(obj, flag) # not sure what flag = "flag" does <- read online, its for 
opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 1)) #, lower = -5.5, upper = 5.5)
sd_report <- sdreport(obj, getJointPrecision = TRUE)

## Start the bootstrap =========================================================
boot_fits <- list()

B <- 10 # number of bootstraps
for(i in 1:B) {
  cat(paste0("Currently on bootstrap number ", i, "of ", B, "\n"))
  set.seed(18012021 + i)
  boot_rows <- sample(1:nrow(data$W), nrow(data$W), replace = TRUE)
  
  boot_data <- list(Y_rec = data$Y_rec[boot_rows, ],
                    Y_grid = data$Y_grid,
                    X = data$X,
                    W = data$W[boot_rows, ],
                    R = data$R[boot_rows, ],
                    A = data$A,
                    trunc_level = data$trunc_level,
                    design_matrix = data$design_matrix) # create data
  
  ## Do bootstrap TMB things
  obj_boot <- MakeADFun(boot_data, parameters, DLL="ascr", silent = FALSE)
  # obj <- normalize(obj, flag) # not sure what flat = "flag" does <- read online, its for 
  opt_boot <- stats::nlminb(obj_boot$par, obj_boot$fn, obj_boot$gr, control = list(trace = 0))
  
  boot_fits[[i]] <- opt_boot
}
opt_boot$par

pars <- sapply(boot_fits, function(fit) fit$par)


## Do bootstrap TMB things
obj_boot <- MakeADFun(boot_data, parameters, DLL="ascr", silent = FALSE)
# obj <- normalize(obj, flag) # not sure what flat = "flag" does <- read online, its for 
opt_boot <- stats::nlminb(obj_boot$par, obj_boot$fn, obj_boot$gr, control = list(trace = 1))#, lower = -5.5, upper = 5.5)
sdreport(obj, getJointPrecision = TRUE)



## =============================================================================


library(Rcpp)
library(matrixStats)

sourceCpp("Scripts/bowhead whales/llkRcpp.cpp")

sourceCpp(code = "
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double test(NumericVector x) {
  NumericVector out = cumprod(x);
  double out2 = rev(out)[0];
  
  return out2;
}
")

test(c(1, 2, 3))
