# --- Step 1: Set parameters ---
theta_true <- 5        # True theta
sigma <- 2             # Known SD of X_i
n <- 1000               # Sample size

# Prior hyperparameters
mu0 <- 0
sigma0 <- 10

# --- Step 2: Simulate data ---
set.seed(123)
x <- rnorm(n, mean = theta_true, sd = sigma)

# Histogram of simulated data
hist(x, 
     col = "lightblue", 
     main = "Histogram of Simulated Data",
     xlab = "x values",
     ylab = "Frequency")

# --- Step 3: Compute posterior parameters ---
x_bar <- mean(x)

sigma_n2 <- 1 / (n / sigma^2 + 1 / sigma0^2)   # Posterior variance
mu_n <- sigma_n2 * (n * x_bar / sigma^2 + mu0 / sigma0^2)  # Posterior mean
sigma_n <- sqrt(sigma_n2)  # Posterior SD

# "Sample mean:"
x_bar 
# "Posterior mean (mu_n):"
mu_n
# "Posterior SD (sigma_n):"
sigma_n

# --- Step 4: Draw posterior samples ---
posterior_samples <- rnorm(1000, mean = mu_n, sd = sigma_n)
# "Posterior mean (numerically):"
mean(posterior_samples)
mu_n
# Histogram of posterior

hist(posterior_samples, breaks = 30, probability = TRUE,
     main = "Posterior Distribution of theta",
     xlab = expression(theta),
     col = "pink")

# Overlay posterior density
theta_vals <- seq(mu_n - 4*sigma_n, mu_n + 4*sigma_n, length.out = 500)

lines(theta_vals, dnorm(theta_vals, mean = mu_n, sd = sigma_n), col = "red", lwd = 2)



# --- Step 5: Highlight 95% posterior interval ---
ci <- qnorm(c(0.05, 0.95), mean = mu_n, sd = sigma_n)
abline(v = ci, col = "blue", lwd = 2, lty = 2)
abline(v = mu_n, col = "red", lwd = 2, lty = 2)

# --- Step 6: Highlight 95% posterior interval (by posterior samples) ---
ci_sample=quantile(posterior_samples,c(0.05, 0.95))
abline(v = ci_sample, col = "orange", lwd = 2, lty = 2)

legend("topright", legend = c("Posterior Density", "Posterior Mean", "95% CI"),
       col = c("red", "red", "blue"), lwd = 2, lty = c(1,2,2))


# "95% posterior interval: 
ci
