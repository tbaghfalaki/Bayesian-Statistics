# Load required packages
library(extraDistr)  # for Half-Cauchy

set.seed(123)

# Observed data
y <- list(
  c(5, 7),
  c(10, 12),
  c(4, 6, 5)
)

n_groups <- length(y)
n_sims <- 10000

# Known within-group variance
sigma2 <- 1

# Sample means and sizes
y_bar <- sapply(y, mean)
n_j <- sapply(y, length)

# Monte Carlo samples for theta_j
theta_samples <- matrix(NA, nrow = n_sims, ncol = n_groups)
tau_samples <- numeric(n_sims)

# Function to sample tau from Half-Cauchy(0,5)
sample_half_cauchy <- function(scale = 5) {
  u <- runif(1)
  scale * tan(pi * (u - 0.5))
}

for (i in 1:n_sims) {
  
  # 1. Sample tau from Half-Cauchy(0,5)
  tau <- abs(sample_half_cauchy(5))
  tau_samples[i] <- tau
  tau2 <- tau^2
  
  # 2. Compute posterior mean and variance of theta_j after integrating out mu
  for (j in 1:n_groups) {
    w_j <- n_j[j]/sigma2
    w_others <- sapply(1:n_groups, function(k) 1/(tau2 + sigma2/n_j[k]))
    post_var <- 1 / (w_j + sum(w_others))
    post_mean <- post_var * (w_j * y_bar[j] + sum(w_others * y_bar))
    theta_samples[i,j] <- rnorm(1, mean = post_mean, sd = sqrt(post_var))
  }
}

# Posterior summaries
theta_means <- apply(theta_samples, 2, mean)
theta_sds <- apply(theta_samples, 2, sd)

# Display results
for (j in 1:n_groups) {
  cat(sprintf("Group %d: Posterior mean = %.2f, Posterior SD = %.2f\n",
              j, theta_means[j], theta_sds[j]))
}

# Optional: plot posterior distributions
par(mfrow=c(1, n_groups))
for (j in 1:n_groups) {
  hist(theta_samples[,j], breaks=50, main=paste("Posterior of theta", j),
       xlab=expression(theta[j]), col="skyblue", border="white")
  abline(v = theta_means[j], col="red", lwd=2)
}
