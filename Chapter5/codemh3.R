set.seed(123)

N <- 10000         # Number of samples
sigma_rw <- 1      # Random walk proposal std dev

# Target density (unnormalized)
target_density <- function(theta) {
  0.2 * exp(-0.5 * theta^2) + 0.8 * exp(-0.5 * (theta - 2)^2)
}

samples <- numeric(N)
samples[1] <- 0  # Initial value

accept_count <- 0

for (t in 2:N) {
  candidate <- samples[t - 1] + rnorm(1, mean = 0, sd = sigma_rw)
  
  # Since proposal symmetric, acceptance ratio simplifies
  alpha <- min(1, target_density(candidate) / target_density(samples[t - 1]))
  
  if (runif(1) < alpha) {
    samples[t] <- candidate
    accept_count <- accept_count + 1
  } else {
    samples[t] <- samples[t - 1]
  }
}

acceptance_rate <- accept_count / (N - 1)
cat("Acceptance rate:", round(acceptance_rate, 4), "\n")

# Trace plot of last 500 samples
plot(samples[(N-499):N], type = "l", col = "steelblue", lwd = 2,
     main = "Trace Plot of Random Walk MH Samples",
     xlab = "Iteration", ylab = expression(theta))

# Histogram with overlays
hist(samples, breaks = 50, probability = TRUE, col = rgb(0, 0, 1, 0.3),
     main = "Histogram of RW MH Samples with Target Density",
     xlab = expression(theta), ylim = c(0, 0.5))

# Approximate normalization constant
norm_const <- integrate(target_density, lower = -Inf, upper = Inf)$value

curve(target_density(x) / norm_const, from = min(samples), to = max(samples),
      add = TRUE, col = "red", lwd = 2)

# Kernel density estimate of samples
lines(density(samples), col = "blue", lwd = 2, lty = 3)

legend("topright", legend = c("Target Density", "Kernel Density"),
       col = c("red", "blue"), lwd = 2, lty = c(1, 3))
