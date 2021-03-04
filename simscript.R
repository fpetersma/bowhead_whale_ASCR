# simulation script test

results <- lapply(c(150:159, 250:259), function(seed) {
  print(paste0("Currently at seed : ", seed))
  
  source("Scripts/Bowhead Whales/simulate_main.R", local = TRUE)
  
  save(file = paste0("temp save with bearings seed = ", seed, ".RData"), 
       list = c("bw"))
  
  return(bw)
})


par(mfrow = c(3, 3))

# If variable sl
plot(density(sapply(results, function(x){x$real[7, 2]})), 
     main = "Distribution of estimates for 'D'")
abline(v = exp(-2), col = "red")

# If single sl
plot(density(sapply(results, function(x){x$real[5, 2]})), 
     main = "Distribution of estimates for 'D'")
abline(v = exp(-2), col = "red")

# Both scenarios
plot(density(sapply(results, function(x){x$real[1, 2]})), 
     main = "Distribution of estimates for 'p'")
abline(v = 0.5, col = "red")

# With bearings ======
plot(density(sapply(results, function(x){x$real[2, 2]})), 
     main = "Distribution of estimates for 'kappa'")
abline(v = 500, col = "red")
# ====================

plot(density(sapply(results, function(x){x$real[3, 2]})), 
     main = "Distribution of estimates for 'beta_r'")
abline(v = 15, col = "red")

plot(density(sapply(results, function(x){x$real[4, 2]})), 
     main = "Distribution of estimates for 'sd_r'")
abline(v = 3, col = "red")

plot(density(sapply(results, function(x){x$real[5, 2]})), 
     main = "Distribution of estimates for 'mu_s'")
abline(v = 140, col = "red")

# if variable sl
plot(density(sapply(results, function(x){x$real[6, 2]})), 
     main = "Distribution of estimates for 'sd_s'")
abline(v = 5, col = "red")



# ON THE LINK =================================

# If variable sl
plot(density(sapply(results, function(x){x$estimates[6, 2]})), 
     main = "Distribution of estimates for 'D'")
abline(v = -2, col = "red")

# If single sl
plot(density(sapply(results, function(x){x$estimates[5, 2]})), 
     main = "Distribution of estimates for 'D'")
abline(v = -2, col = "red")

# Both scenarios
plot(density(sapply(results, function(x){x$estimates[1, 2]})), 
     main = "Distribution of estimates for 'p'")
abline(v = log(0.5 / (1 - 0.5)), col = "red")

plot(density(sapply(results, function(x){x$estimates[2, 2]})), 
     main = "Distribution of estimates for 'beta_r'")
abline(v = log(15), col = "red")

plot(density(sapply(results, function(x){x$estimates[3, 2]})), 
     main = "Distribution of estimates for 'sd_r'")
abline(v = log(3), col = "red")

plot(density(sapply(results, function(x){x$estimates[4, 2]})), 
     main = "Distribution of estimates for 'mu_s'")
abline(v = log(140), col = "red")

# if variable sl
plot(density(sapply(results, function(x){x$estimates[5, 2]})), 
     main = "Distribution of estimates for 'sd_s'")
abline(v = log(5), col = "red")
