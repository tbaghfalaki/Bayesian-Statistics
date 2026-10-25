set.seed(123)  # for reproducibility

# Parameters
mu <- 2
sigma <- 1
n_samples <- 10000

# Rejection sampling constants
k <- (2 / sqrt(2 * pi)) * exp(0.5)  # approx 1.315

# Function to compute acceptance probability
accept_prob <- function(x) {
  return( exp(-(x - 1)^2 / 2) / k )
}

# Rejection sampling from half-normal (|Z|)
sample_half_normal <- function(n) {
  samples <- numeric(n)
  i <- 1
  while (i <= n) {
    # Sample from exponential(1)
    X <- rexp(1, rate = 1)
    U <- runif(1)
    if (U <= accept_prob(X)) {
      samples[i] <- X
      i <- i + 1
    }
  }
  return(samples)
}

# Generate half-normal samples
X_samples <- sample_half_normal(n_samples)


# Transform to Y ~ N(mu, sigma^2)
Y_samples <- mu + sigma * Z_samples

# Plot histogram of Y_samples with theoretical normal density
hist(Y_samples, breaks = 50, probability = TRUE,
     main = "Histogram of Samples from HN(μ, σ²) via Rejection Sampling",
     xlab = "Y", col = "lightblue", border = "white")

# Add theoretical normal density curve
curve(2*dnorm(x, mean = mu, sd = sigma), col = "red", lwd = 2, add = TRUE)
legend("topright", legend = c("Sample histogram", "True Half Normal density"),
       col = c("lightblue", "red"), lwd = c(10, 2), bty = "n")