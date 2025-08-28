set.seed(42)

# Parameters
N <- 10000          # Number of samples
sigma_proposal <- 1 # Proposal standard deviation

# Target density: standard Cauchy (unnormalized)
cauchy_density <- function(x) {
  1 / (pi * (1 + x^2))
}

# Proposal density: Normal centered at current state
proposal_density <- function(x, mean) {
  dnorm(x, mean = mean, sd = sigma_proposal)
}

# Initialize chain
samples <- numeric(N)
samples[1] <- 0  # initial value

accept_count <- 0

for (t in 2:N) {
  candidate <- rnorm(1, mean = samples[t - 1], sd = sigma_proposal)
  
  # Symmetric proposal => acceptance ratio = target(candidate)/target(current)
  alpha <- min(1, cauchy_density(candidate) / cauchy_density(samples[t - 1]))
  
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
plot(samples[(N-499):N], type = "l", lwd = 2, col = "steelblue",
     ylab = "Sample value", xlab = "Iteration",
     main = "Trace Plot of MH Samples from Cauchy(0,1) with Normal(0,1) Proposal")


# Plot histogram and densities
hist(samples, breaks = 50, probability = TRUE, col = rgb(0, 0, 1, 0.3),
     main = "MH Sampling from Cauchy(0,1) with Normal(0,1) Proposal",
     xlab = expression(paste(theta)), ylim = c(0, 0.5))

# Overlay true target density (solid red)
curve(dcauchy(x), add = TRUE, col = "red", lwd = 2)

# Overlay proposal density (dashed green)
# Proposal density is Normal(0,1), but the proposal centers vary at each step,
# so for visualization, use Normal(0,1) fixed:
curve(dnorm(x, mean = 0, sd = sigma_proposal), add = TRUE, col = "darkgreen",
      lwd = 2, lty = 2)

# Overlay kernel density estimate of samples (dotted blue)
lines(density(samples), col = "blue", lwd = 2, lty = 3)

legend("topright", legend = c("True Cauchy Density", "Proposal Density", "Kernel Density"),
       col = c("red", "darkgreen", "blue"), lwd = 2, lty = c(1, 2, 3))
