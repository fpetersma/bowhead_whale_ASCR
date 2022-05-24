# dat <- data_obj; par <- parameters; rm(list = setdiff(ls(), c("dat", "par")));
## Create R wrapper
nllR <- function(par, dat) {
  ## Print the parameters
  # print(round(par, 5))
  # cat("current pars: ", par, "\n", sep = " ")
  
  ## Write table for parameter history
  # write.table(matrix(par, nrow = 1, byrow = TRUE),  
  #             file = paste0("parameter_history_Rcpp_", dat$f_density[3], ".csv"),
  #             append = T, 
  #             sep = ',', 
  #             row.names = FALSE, 
  #             col.names = FALSE)
  
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
  gam_fit <- dat$gam_fit
  
  ## Derive density based on GAM object
  gam_fit$coefficients <- pars$density_pars         # match order of parameters!
  eta <- as.vector(predict(gam_fit, data = model.matrix(gam_fit)))
  dat$D <- exp(eta) # linear in the exponent to ensure positive density
  
  if (dat$USE_BEARINGS == 1) {
    ## Get the besselI function returns
    dat$bessel <- besselI(exp(pars$log_kappa), 0)
  }
  if (dat$USE_BEARINGS == 2) {
    ## Get the besselI function returns
    dat$bessel_low <- besselI(exp(pars$log_kappa_low), 0)
    dat$bessel_high <- besselI(exp(pars$log_kappa_low) + 
                                 exp(pars$log_kappa_high), 0)
  }
  
  ## ===========================================================================
  ## Run below to compare simpson and midpoint integration
  ## ---------------------------------------------------------------------------
  ## Simpson's rule
  
  # sourceCpp("Scripts/bowhead whales/nllRcppFlexBearSimpson.cpp")
  # integration <- "simpson"
  # if (integration == "simpson") {
  #   dat$n_sl <- 30 + 1
  #   dat$S <- seq(from = 80, to = 240, length.out = dat$n_sl)  
  #   dat$A_s <- dat$S[2] - dat$S[1]
  #   S <- dat$S
  # }
  # if (!dat$FIXED_SL) {
  #   ## Probability of calls (potentially truncated normal?)
  #   dat$S_probs_log <- dnorm(S, exp(pars$log_mu_s), exp(pars$log_sd_s), TRUE)
  #   dat$S_probs <- exp(dat$S_probs_log)
  # }
  # uncondNllRcppVarSL(dat, par)
  
  ## Midpoint rule
  
  # sourceCpp("Scripts/bowhead whales/nllRcppFlexBear.cpp")
  # integration <- "midpoint"
  # if (integration == "midpoint") {
  #   dat$n_sl <- 30
  #   dat$S <- seq(from = 100, to = 220, length.out = dat$n_sl + 1)
  #   dat$A_s <- diff(dat$S)[1]
  #   dat$S <- dat$S[-length(dat$S)] + dat$A_s / 2
  #   S <- dat$S
  # }
  # if (!dat$FIXED_SL) {
  #   ## Probability of calls (potentially truncated normal?)
  #   dat$S_probs_log <- dnorm(S, exp(pars$log_mu_s), exp(pars$log_sd_s), TRUE)
  #   dat$S_probs <- exp(dat$S_probs_log)
  # }
  # uncondNllRcpp(dat, par)
  ## ===========================================================================
  
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
  ## Parallelise on the calls if SINGLE_SL == FALSE
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
      # return(uncondNllRcppVarSL(dat_i, par))
    })
  } else {
    cores <- 6
    
    ## Initiate parallel process leaving one core free
    no_cores <- cores #detectCores() - 1
    cl <- makeCluster(no_cores)
    
    ## Export relevant data to all clusters
    clusterExport(cl, varlist = c("dat", "par"), envir = environment())
    
    ## Start the parallelisation 
    nll_cond_call <- parLapply(cl, as.list(1:dat$n_call), function(i) {
      ## Extract the data only for call i
      dat_i <- dat
      dat_i$Y_rec <- dat_i$Y_rec[i, ]
      dat_i$W <- dat_i$W[i, ]
      dat_i$R <- dat_i$R[i, ]
      
      # Rcpp::sourceCpp("Scripts/bowhead whales/nllRcppFlexBearSimpson.cpp")
      return(ascrRcpp::singleNllRcpp(dat_i, par))
      # return(singleNllRcppVarSL(dat_i, par))
    })
    
    stopCluster(cl)
  }
  
  nll_cond <- sum(unlist(nll_cond_call))
  
  ## ===========================================================================
  
  ## ===========================================================================
  ## sum it all together
  ## ---------------------------------------------------------------------------
  nll <- nll_uncond + nll_cond
  
  ## ===========================================================================
  
  # cat("Current nll is: ", nll, "\n", sep = "")
  
  ## Return the negative log likelihood
  return(nll)
}
