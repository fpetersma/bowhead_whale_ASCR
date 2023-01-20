# dat <- data_obj; par <- parameters; rm(list = setdiff(ls(), c("dat", "par")));
## Create R wrapper
nllR <- function(par, dat) {
  ## Print the parameters
  # print(round(par, 5))
  # cat("current pars: ", par, "\n", sep = " ")
# 
#   ## Write table for parameter history [[don't use for bootstrapping]]
#   write.table(matrix(par, nrow = 1, byrow = TRUE),
#               file = paste0("parameter_history_Rcpp_", dat$f_density[3], ".csv"),
#               append = T,
#               sep = ',',
#               row.names = FALSE,
#               col.names = FALSE)
#   
  ## ===========================================================================
  ## Extract parameters and data
  ## ---------------------------------------------------------------------------
  
  pars <- as.list(par[!(names(par) %in% dat$dens_par_names)])
  pars[["density_pars"]] <- par[names(par) %in% dat$dens_par_names]
  
  S <- dat$S
  # design_matrix <- data$design_matrix
  # 
  ## ===========================================================================
  
  
  ## ===========================================================================
  ## DERIVE INFORMATION THAT ONLY DEPENDS ON CURRENT PARAMETER VALUES
  ## ---------------------------------------------------------------------------
  gam_fit_weird <- dat$gam_fit
  
  ## Derive density based on GAM object
  gam_fit_weird$coefficients <- pars$density_pars         # match order of parameters!
  eta <- as.vector(mgcv::predict.gam(gam_fit_weird, data = model.matrix(gam_fit_weird)))
  dat$D <- exp(eta) # linear in the exponent to ensure positive density
  dat$gam_fit <- gam_fit_weird
  
  if (dat$USE_BEARINGS == 1) {
    ## Get the besselI function returns
    dat$bessel <- besselI(exp(pars$log_kappa), 0)
  }
  if (dat$USE_BEARINGS == 2) {
    ## Get the besselI function returns (limit it to 1e100)
    dat$bessel_low <- min(besselI(exp(pars$log_kappa_low), 0), 1e100)
    dat$bessel_high <- min(besselI(exp(pars$log_kappa_low) + 
                                   exp(pars$log_kappa_high), 0), 1e100)
    
  }

  
  if (!dat$FIXED_SL) {
    ## Probability of calls (potentially truncated normal?)
    dat$S_probs_log <- dnorm(S, exp(pars$log_mu_s), exp(pars$log_sd_s), TRUE)
    dat$S_probs <- exp(dat$S_probs_log)
  }
  
  # sourceCpp("Scripts/bowhead whales/nllRcppFlexBear.cpp")
  
  ## ===========================================================================
  ## COMPLETE THE UNCONDITIONAL LOG LIKELIHOOD PART
  ## ---------------------------------------------------------------------------
  if (dat$FIXED_SL) {
    nll_uncond <- ascrRcpp::uncondNllRcppFixedSL(dat, par)
  } else {
    nll_uncond <- ascrRcpp::uncondNllRcpp(dat, par)
  }
  
  ## ===========================================================================
  
  ## ===========================================================================
  ## Parallelise on the calls if FIXED_SL == FALSE
  ## ---------------------------------------------------------------------------
  
  if (dat$FIXED_SL) {
    ## Start the parallelisation 
    nll_cond_call <- lapply(as.list(1:dat$n_call), function(i) {
      ## Extract the data only for call i
      dat_i <- dat
      dat_i$Y_rec <- dat_i$Y_rec[i, ]
      dat_i$W <- dat_i$W[i, ]
      dat_i$R <- dat_i$R[i, ]
      
      return(ascrRcpp::singleNllRcppFixedSL(dat_i, par))
    })
  } else {
    # cores <- 4
    # 
    # ## Initiate parallel process leaving one core free
    # no_cores <- cores #detectCores() - 1
    # cl <- parallel::makeCluster(no_cores)
    # 
    # ## Export relevant data to all clusters
    # parallel::clusterExport(cl, varlist = c("dat", "par"), envir = environment())
    # 
    # ## Start the parallelisation
    # nll_cond_call <- parallel::parLapply(1:dat$n_call, function(i) {
    #   ## Extract the data only for call i
    #   dat_i <- dat
    #   dat_i$Y_rec <- dat_i$Y_rec[i, ]
    #   dat_i$W <- dat_i$W[i, ]
    #   dat_i$R <- dat_i$R[i, ]
    # 
    #   return(ascrRcpp::singleNllRcpp(dat_i, par))
    # }, cl = cl); parallel::stopCluster(cl)

    # Non parallel version
    nll_cond_call <- lapply(1:dat$n_call, function(i) {
    # nll_cond_call <- lapply(1:50, function(i) {
      ## Extract the data only for call i
      dat_i <- dat
      dat_i$Y_rec <- dat_i$Y_rec[i, ]
      dat_i$W <- dat_i$W[i, ]
      dat_i$R <- dat_i$R[i, ]

      return(ascrRcpp::singleNllRcpp(dat_i, par))
    })
  }
  
  nll_cond <- sum(unlist(nll_cond_call))
  
  ## ===========================================================================
  
  ## ===========================================================================
  ## sum it all together
  ## ---------------------------------------------------------------------------
  nll <- nll_uncond + nll_cond
  
  ## ===========================================================================
  
  # cat("Current nll is:", nll, "and current pars are:\n", paste0(names(par), ": ", par, "\n"), "\n", sep = " ")
  
  ## Return the negative log likelihood
  return(nll)
  
  ## Return the log likelihood
  # return(-nll)
}
