# A function that calculates the density
.densityGAM <- function(x, gam_fit, par) {
  # Description:
  #   Takes in a matrix of (x, y)-coordinates and calculates the density at 
  #   every one of them. So far only working for matrices of M x 2 format, for
  #   M is mesh size. 
  
  # Inputs:
  #   x - x & y coordinates of the grid and any additional spatial covariates
  #   f - the specified relation between D and the spatial covariates
  #   param - a vector of parameters corresponding to the provided density
  #     relation.
  
  # Outputs:
  #   D - a vector of densities per unit area
  
  # Match order to parameters
  gam_fit$coefficients <- par
  
  eta <- predict(gam_fit, data = model.matrix(gam_fit))
  
  D <- exp(eta) # linear in the exponent to ensure positive density
  
  if (any(D < 0)) {
    stop("Negative estimates for density!")
  }
  
  return(as.vector(D))
}

.gSound <- function(level, par, type = "simple", ...) {
  # Description:
  #   
  
  # Inputs:
  #   
  # Outputs:
  args <- list(...)

  if (type == "simple") {
    g0 <- par["g0"]
    g <- g0 * (1 - pnorm((args$trunc_level - level) / args$sd_r)) # finds sd_r and trunc_level in ...
  } else {
    U <- par["U"]
    B <- par["B"]
    Q <- par["Q"]
    L <- 0
    
    if (type == "probit") {
      if (U < 0 | U > 1) {stop("The upper limit U should be between 0 and 1!")}
      # if (B <= 0) {stop("The inflection point parameter B should be positive!")}
      if (Q <= 0) {stop("The growth parameter Q should be larger than 1!")}
      
      g <- U * pnorm(level, B, Q)
    }
    if (type == "logit") {
      # Input checks
      if (U < 0 | U > 1) {stop("The upper limit U should be between 0 and 1!")}
      if (B <= 0) {stop("The inflection point parameter B should be positive!")}
      if (Q <= 0) {stop("The growth parameter Q should be larger than 1!")}
      
      g <- U * plogis(level, B, Q)
    }
    if (type == "janoschek") {
      # Input checks
      if (U < 0 | U > 1) {stop("The upper limit U should be between 0 and 1!")}
      if (B <= 0) {stop("The inflection point parameter B should be positive!")}
      if (Q <= 1) {stop("The growth parameter Q should be larger than 1!")}
      
      g <- U - (U - L) * exp(-(B / 1000) * (level ^ Q)) # B divided by 1000 to improve convergence
    }
  }
  
  return(g)
}

.gJano <- function(snr, par) {
  # Description:
  #   
  
  # Inputs:
  #   
  # Outputs:
  
  U <- par["U"]
  B <- par["B"]
  Q <- par["Q"]
  L <- 0
  
  # Input checks
  if (U < 0 | U > 1) {stop("The upper limit U should be between 0 and 1!")}
  if (B <= 0) {stop("The inflection point parameter B should be positive!")}
  if (Q <= 1) {stop("The growth parameter Q should be larger than 1!")}
  
  
  g <- U - (U - L) * exp(-(B / 1000) * (snr ^ Q)) # B divided by 1000 to improve convergence

  return(g)
}

.gHN <- function(distances, par) {
  g0 <- par["g0"]
  sigma <- par["sigma"]
  
  if (g0 <= 0) {stop("The g0 parameter B should be positive!")}
  if (sigma <= 0) {stop("The range parameter sigma should be positive!")}
  
  g <- g0 * exp(- distances ^ 2 / (2 * sigma ^ 2))
  
  return(g)
}

# Using column operations (should be faster)
.detected <- function(probs, min_no_detections) {
  
  # Description:
  
  # Inputs:
  
  # Outputs:
  
  
  probs <- t(probs)
  
  k <- nrow(probs)
  
  # Create all combinations of detections and non-detections
  dummies <- t(expand.grid(rep(list(1:0), k)))
  
  # Only keep row with less than 'min_no_detections' detections
  dummies <- dummies[, Rfast::colsums(dummies) < min_no_detections]
  
  # Calculate the probabilities of detection animals less than the min_no_detections
  temp <- apply(probs, 2, function(p1) {
    p0 <- 1 - p1
    out <- matrix(rep(p0, ncol(dummies)), nrow = nrow(dummies), byrow = FALSE)
    p1_m <- matrix(rep(p1, ncol(dummies)), nrow = nrow(dummies), byrow = FALSE)
    out[as.logical(dummies)] <- p1_m[as.logical(dummies)]
    # Take the product of all probabilities
    return(sum(Rfast::colprods(out)))
    
  })
  # Calculate probability that it was detected at least twice
  p. <- 1 - temp
  
  # When probabiltiies of detection are really small, the probability of not
  # seen is 1, such that 1 - 1 - P(seen once) < 0. This is of course impossible,
  # a probability can never be negative, so I now simply set everything < 0
  # equal to zero. Not sure if this is a very nice way to handle this.
  p.[p. < 0] <- 0
  
  return(p.)
}

.totalN <- function(D, A, covariance_matrix, design, n) {
  # Description:
  # Inputs:
  #   D - vector of density corresponding to A
  #   A - area of every grid cell corresponding to D
  #   est - estimates corresponding to the densit specification
  #   n - number of calls detected at least twice
  # Input checks!
  # D <- .density_GAM(design, f, est)
  N <- sum(A * D)
  
  se <- NULL
  upper <- NULL
  lower <- NULL
  
  try({
    # Calculate the variance of N based on Millar (2011) or Murray (2013). 
    # g(x) is defined as sum(D) [see notes]
    density_derivatives <- unname(colSums(D * design))
    G <- c(rep(0, length = ncol(covariance_matrix) - length(density_derivatives)), 
           density_derivatives)
    v <- t(G) %*% covariance_matrix %*% G
    v <- mean(A) ^ 2 * as.vector(v) - N # mean(A) is almost surely wrong, and the - N is based on eq (4) in Efford et al 2013.
    
    # Assume normal distribution of N 
    se <- sqrt(v)  
    
    # Since we know the number of detected animals for sure (= n), we don't need
    # to include that part of the expectation in the variance. Basically, there is
    # only uncertainty in N - n. Use Murray (2013) method for lognormal interval,
    # in line with Rexstad and Burnham (1991).
    C <- exp(1.96 * sqrt(log(1 + (v / (N - n)) ^ 2)))
    
    lower <- n + (N - n) / C
    upper <- n + (N - n) * C
  })
  return(c(estimate = round(N, 6), 
           se = round(se, 6), 
           lower = round(lower, 6), 
           upper = round(upper, 6)))
}
