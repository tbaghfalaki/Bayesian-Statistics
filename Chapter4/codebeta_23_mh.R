# Target distribution parameters
a <- 2
b <- 3
N <- 5000  # Number of samples

# Initialize vectors
X <- numeric(N)
alpha <- numeric(N)

# Initial value for the chain (sampled from Uniform(0,1))
X[1] <- runif(1)

# Metropolis-Hastings sampling
for (n in 2:N) {
  proposal <- runif(1)  # Propose a new value from Uniform(0,1)
  
  # Compute acceptance probability
  alpha[n] <- min(1, dbeta(proposal, a, b) / dbeta(X[n - 1], a, b))
  
  # Accept/reject the proposal
  if (runif(1) < alpha[n]) {
    X[n] <- proposal  # Accept the proposal
  } else {
    X[n] <- X[n - 1]  # Reject: stay at the current state
  }
}

par(mfrow=c(2,1))
# Trace plot of the last 500 samples to assess mixing
plot(X[(N - 499):N], type = "l", lwd = 2, col = "steelblue",
     ylab = "Sample value", xlab = "Iteration",
     main = "Trace Plot of MH Samples from Beta(2,3)")

# Histogram with true density overlaid
hist(X, breaks = 40, probability = TRUE, col = "lightblue",
     main = "Histogram of MH Samples with Beta(2,3) Density",
     xlab = expression(theta))
curve(dbeta(x, a, b), add = TRUE, col = "red", lwd = 2)

# Report the mean acceptance rate
cat("Mean acceptance rate:", round(mean(alpha[-1]), 3), "\n")
