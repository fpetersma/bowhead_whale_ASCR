bufferWidth <- function(p,
                        trunc_level,
                        g0, 
                        beta_r, 
                        sd_r, 
                        mu_s, 
                        sd_s = NULL, 
                        noise = NULL,
                        SINGLE_SL = TRUE, 
                        WITH_NOISE = FALSE) {
  # Get the correct quantile for the normal cdf, given p
  quant <- qnorm(1 - p)
  
  # Derive the distance at which that p is achieved
  distance <- 10 ^ ((quant * sd_r - trunc_level + mu_s) / beta_r)
  
  return(distance)
  # Assuming all detectors are at this distance, what is p.? 
  # RECALL: we are only interested in calls that were detected at least twice!)
  
}

bufferWidth(p = 0.001,
            trunc_level = 96,
            g0 = 0.62,
            beta_r = 17.36,
            sd_r = 4.85,
            mu_s = 166.3)
