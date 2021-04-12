# Run for test runs ============================================================
# rm(list = setdiff(ls(), c("dat", "par_fix", "par_flex")))
# TRACE <- FALSE
# USE_BEARINGS <- TRUE
# ==============================================================================

llkSurface <- function(par_fix, par_flex, dat, USE_BEARINGS) {
  
  # ============================================================================
  # DO ALL THE PREWORK BASED ON bwASCR() 
  # ============================================================================
  
  
  ## Load libraries and perform input checks -----------------------------------
  library(dplyr)
  library(mgcv)
  library(matrixStats)
  library(circular)
  library(raster)
  
  # Create a matrix of distances
  distances <- pointDistance(p1 = dat$cov_density[, c("long", "lat")], 
                             p2 = dat$detectors[, c("long", "lat")], 
                             lonlat = TRUE)
  dat[["distances"]] <- distances  # Add matrix of distances to dat
  
  # Extract variables from f_density, remove 'D' and add '(Intercept)'
  vars <- all.vars(dat$f_density)
  vars <- c("(Intercept)", vars[-1])
  
  # Check whether density is constant
  if (length(vars) == 1) {
    CONSTANT_DENSITY <- FALSE
  } else {
    CONSTANT_DENSITY <- FALSE
  }
  dat["CONSTANT_DENSITY"] <- CONSTANT_DENSITY
  dat[["USE_RL"]] <- TRUE
  dat[["TRACE"]] <- FALSE
  dat[["USE_BEARINGS"]] <- USE_BEARINGS
  
  # Not sure if scaling is necessary
  cov_density_scaled <- subset(dat$cov_density, select = -area)
  cov_density_scaled[, -c(1, 2)] <- scale(subset(cov_density_scaled, 
                                                 select = -c(long, lat)))
  
  # Turn the bearings to radians, as the llk calculation uses radians 
  if (USE_BEARINGS) {
    bearings_deg <- circular(dat$bearings,
                             units = "degrees",
                             template = "geographics") # 'geographics' means clockwise rotation with 0 at north
    bearings_rad <- conversion.circular(bearings_deg)
    dat[["bearings_rad"]] <- bearings_rad # Add bearings as circular in radians to dat
    
    # Get degrees bearings from detectors to all grid points 
    grid_bearings <- t(apply(dat$cov_density[, c("long", "lat")], 1, function(grid_point) {
      bear <- geosphere::bearing(p1 = as.matrix(dat$detectors[, c("long", "lat")]),
                                 p2 = grid_point)
    }))
    # First convert to 'circular' 
    grid_bearings <- circular(grid_bearings, template = "geographics", 
                              modulo = "2pi", units = "degrees")
    # Then convert to radians
    grid_bearings <- conversion.circular(grid_bearings, units = "radians")
    dat[["grid_bearings"]] <- grid_bearings
  }
  
  if (length(par_flex) == 1) {
    
  } else if (length(par_flex) == 2) {
    llk_values <- sapply(par_flex[[1]], function(par_flex_1) {
      flex2 <- sapply(par_flex[[2]], function(par_flex_2) {
        names(par_flex_1) <- names(par_flex)[1]
        names(par_flex_2) <- names(par_flex)[2]
        
        # Turn par into a named vector with the correct names (optim() requires a vector)
        par <- c(par_flex_1, par_flex_2, unlist(par_fix))
        names(par) <- gsub(".*\\.", "", names(par))
        
        if (dat$SINGLE_SL) {
          par <- par[!names(par) %in% c("sd_s")] # USE ONLY FOR SINGLE SOURCE LEVEL
        }
        
        
        # Create a gam object (not necessary if CONSTANT_DENSITY == FALSE?)
        gam_fit <- gam(dat$f_density, 
                       data = cbind(D = 0, cov_density_scaled))
        if (length(gam_fit$smooth) != 0) {
          k <- gam_fit$smooth[[1]]$bs.dim
          gam_par <- rep(0, k - 1)
          dat[["smooth_terms"]] <- paste0(gam_fit$smooth[[1]]$label, ".", 1:(k - 1))
          names(gam_par) <- dat$smooth_terms
          
          # par <- c(par, gam_par)
          par <- c(par, gam_fit$coefficients[names(gam_fit$coefficients) != "(Intercept)"])
        } else {dat[["smooth_terms"]] <- NULL}
        
        dat[["design"]] <- model.matrix(gam_fit)
        dat[["gam_fit"]] <- gam_fit
        
        return(.llkParallelSmooth(par, dat))
      })
      return(flex2)
    })
    colnames(llk_values) <- round(par_flex[[1]], 4)
    rownames(llk_values) <- round(par_flex[[2]], 4)
  } else {
    stop("number of columns of par_flex has to be either 1 or 2!")
  }
  
  return(llk_values)

  
  
  
  
  
  
  
  # ============================================================================
  # Derive the likelihood values 
  # ============================================================================
  
  
  
}