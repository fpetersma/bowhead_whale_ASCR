## Load the .RData files and give them new names

par_cond <- t(sapply(bw_list_cond, function(x) x$result$par))
par_full <- t(sapply(bw_list_full, function(x) x$result$par))

logit <- function(x, INV = FALSE) {
  if (INV) {
    return(exp(x) / (1 + exp(x)))
  } else {
    return(log(x / (1 - x)))
  }
}

# Compare density from conditional and full likelihood
plot(density(sapply(bw_list_cond, function(x) x$constant_density)), type = "n")
abline(v = exp(-2), col = "red")
lines(density(sapply(bw_list_cond, function(x) x$constant_density)), lty = 2, col = "green")
lines(density(exp(sapply(bw_list_full, function(x) x$result$par[["(Intercept)"]]))), 
      lty = 3, col = "blue")
# They are identical!

plot(density(logit(par_cond[, "g0"], TRUE)), type = 'n')
abline(v = 0.7, col = 'red')
lines(density(logit(par_cond[, "g0"], TRUE)), lty = 2, col = 'green')
lines(density(logit(par_full[, "g0"], TRUE)), lty = 3, col = 'blue')

plot(density(exp(par_cond[, "sigma"])), type = 'n')
abline(v = 15000, col = 'red')
lines(density(exp(par_cond[, "sigma"])), lty = 2, col = 'green')
lines(density(exp(par_full[, "sigma"])), lty = 3, col = 'blue')

plot(density(exp(par_cond[, "kappa"])), type = 'n')
abline(v = 10, col = 'red')
lines(density(exp(par_cond[, "kappa"])), lty = 2, col = 'green')
lines(density(exp(par_full[, "kappa"])), lty = 3, col = 'blue')

### Load some new data, for example the cond and full likelihood version of janoschek
## NOTE: I know here that my density estimate for the cond likelihood is wrong!

bw_list_cond <- bw_list
bw_list_full <- bw_list
rm(bw_list)

par_cond <- t(sapply(bw_list_cond, function(x) x$result$par))
par_full <- t(sapply(bw_list_full, function(x) x$result$par))

# Plot the intercept
plot(density(exp(par_full[, "(Intercept)"])))
abline(v = exp(-1.5), col = "red")

# Plot U
plot(density(logit(par_cond[, "U"], TRUE)), type = 'n')
abline(v = 0.8, col = 'red')
lines(density(logit(par_cond[, "U"], TRUE)), lty = 2, col = 'green')
lines(density(logit(par_full[, "U"], TRUE)), lty = 3, col = 'blue')
# Identical!

# Plot B
plot(density(exp(par_cond[, "B"])), type = 'n')
abline(v = 2, col = 'red')
lines(density(exp(par_cond[, "B"])), lty = 2, col = 'green')
lines(density(exp(par_full[, "B"])), lty = 3, col = 'blue')
# Identical!

# Plot Q
plot(density(exp(par_cond[, "Q"])), type = 'n')
abline(v = 3, col = 'red')
lines(density(exp(par_cond[, "Q"])), lty = 2, col = 'green')
lines(density(exp(par_full[, "Q"])), lty = 3, col = 'blue')
# Identical!

# Plot kappa
plot(density(exp(par_cond[, "kappa"])), type = 'n')
abline(v = 10, col = 'red')
lines(density(exp(par_cond[, "kappa"])), lty = 2, col = 'green')
lines(density(exp(par_full[, "kappa"])), lty = 3, col = 'blue')
# Identical!

# Plot mu_s
plot(density(exp(par_cond[, "mu_s"])), type = 'n')
abline(v = 100, col = 'red')
lines(density(exp(par_cond[, "mu_s"])), lty = 2, col = 'green')
lines(density(exp(par_full[, "mu_s"])), lty = 3, col = 'blue')
# Identical!

# Plot sd_s
plot(density(exp(par_cond[, "sd_s"])), type = 'n')
abline(v = 5, col = 'red')
lines(density(exp(par_cond[, "sd_s"])), lty = 2, col = 'green')
lines(density(exp(par_full[, "sd_s"])), lty = 3, col = 'blue')
# Identical!

# Plot beta_r
plot(density(exp(par_cond[, "beta_r"])), type = 'n')
abline(v = 15, col = 'red')
lines(density(exp(par_cond[, "beta_r"])), lty = 2, col = 'green')
lines(density(exp(par_full[, "beta_r"])), lty = 3, col = 'blue')
# Identical!

# Plot sd_r
plot(density(exp(par_cond[, "sd_r"])), type = 'n')
abline(v = 1, col = 'red')
lines(density(exp(par_cond[, "sd_r"])), lty = 2, col = 'green')
lines(density(exp(par_full[, "sd_r"])), lty = 3, col = 'blue')
# Identical!

