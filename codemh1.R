set.seed(123)

N <- 10000  # Number of samples

# Target density (unnormalized mixture)
target_density <- function(theta) {
  0.2 * exp(-0.5 * theta^2) + 0.8 * exp(-0.5 * (theta - 2)^2)
}

# Proposal density: Normal(2, 2^2)
proposal_density <- function(theta) {
  dnorm(theta, mean = 2, sd = 2)
}

# Proposal sampler
proposal_sample <- function() {
  rnorm(1, mean = 2, sd = 2)
}

samples <- numeric(N)
samples[1] <- 0  # initial value

accept_count <- 0

for (t in 2:N) {
  candidate <- proposal_sample()
  
  numerator <- target_density(candidate) * proposal_density(samples[t - 1])
  denominator <- target_density(samples[t - 1]) * proposal_density(candidate)
  alpha <- min(1, ifelse(denominator == 0, 1, numerator / denominator))
  
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
     main = "Trace Plot of MH Samples",
     xlab = "Iteration", ylab = expression(theta))

# Histogram with overlay of target and proposal densities
hist(samples, breaks = 50, probability = TRUE, col = rgb(0, 0, 1, 0.3),
     main = "Histogram of MH Samples with Target and Proposal Densities",
     xlab = expression(theta), ylim = c(0, 0.5))

# Approximate normalization constant for target density
norm_const <- integrate(target_density, lower = -Inf, upper = Inf)$value

curve(target_density(x) / norm_const, from = min(samples), to = max(samples),
      add = TRUE, col = "red", lwd = 2)

curve(proposal_density(x), from = min(samples), to = max(samples),
      add = TRUE, col = "darkgreen", lwd = 2, lty = 2)

lines(density(samples), col = "blue", lwd = 2, lty = 3)

legend("topright", legend = c("Target Density", "Proposal Density", "Kernel Density"),
       col = c("red", "darkgreen", "blue"), lwd = 2, lty = c(1, 2, 3))
