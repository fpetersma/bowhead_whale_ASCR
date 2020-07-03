### what are the directions of the biases in full likelihood snr det fucntion
#### model?

## Load the .RData files and give them new names

par_full <- t(sapply(bw_list, function(x) x$result$par))

logit <- function(x, INV = FALSE) {
  if (INV) {
    return(exp(x) / (1 + exp(x)))
  } else {
    return(log(x / (1 - x)))
  }
}

# Plot the intercept
plot(density(exp(par_full[, "(Intercept)"])))
abline(v = exp(-1.5), col = "red")
################################## density is slightly overestimated

# Plot U
plot(density(logit(par_full[, "U"], TRUE)), lty = 3, col = 'blue')
abline(v = 0.8, col = 'red')
########################## U is very slightly underestimated


# Plot B
plot(density(exp(par_full[, "B"])), lty = 3, col = 'blue')
abline(v = 2, col = 'red')
######################## B is a little underestimated!

# Plot Q
plot(density(exp(par_full[, "Q"])), lty = 3, col = 'blue')
abline(v = 3, col = 'red')
######################## Very slightly underestimation of Q

# Plot kappa
plot(density(exp(par_full[, "kappa"])), lty = 3, col = 'blue')
abline(v = 10, col = 'red')
####  Underestimated kappa, look into this!!!

# Plot mu_s
plot(density(exp(par_full[, "mu_s"])), lty = 3, col = 'blue')
abline(v = 100, col = 'red')
########################3 Very slightly underestimated mu_s

# Plot sd_s
plot(density(exp(par_full[, "sd_s"])), lty = 3, col = 'blue')
abline(v = 5, col = 'red')
######################### Very slightly overestimated sd_s

# Plot beta_r
plot(density(exp(par_full[, "beta_r"])), lty = 3, col = 'blue')
abline(v = 15, col = 'red')
################################# Slightly underestimated beta_r

# Plot sd_r
plot(density(exp(par_full[, "sd_r"])), lty = 3, col = 'blue')
abline(v = 1, col = 'red')
########################### Overestimated sd_r

###################3 What happens if I do a run where I fix kappa, sd_r and sd_s?
# Start:
#   
#   U           B           Q      beta_r        mu_s (Intercept) 
# 1.3862944   0.6931472   1.0986123   2.7080502   2.9957323  -1.5000000 
# 
# End:
# U           B           Q      beta_r        mu_s (Intercept) 
# 1.346857    1.358784    1.972704    2.922408    2.995732   -1.023466

# I need to check this, cuz this does not feel good. But not now.