library(mvtnorm)

# ---- Target Distribution ----
mean_target <- c(0, 0)
cov_target <- matrix(c(1, -0.6, -0.6, 1), nrow = 2)

# ---- Proposal: Random Walk with Diagonal Covariance ----
proposal_sd <- 0.4  # standard deviation for each coordinate
cov_proposal <- diag(proposal_sd^2, 2)  # diagonal proposal covariance

# ---- MCMC Settings ----
N <- 5000
samples <- matrix(0, nrow = N, ncol = 2)
samples[1, ] <- c(0, 0)
accept_count <- 0

# ---- Random Walk Metropolis-Hastings Algorithm ----
for (t in 2:N) {
  # Propose from N(current, proposal_cov)
  candidate <- rmvnorm(1, mean = samples[t - 1, ], sigma = cov_proposal)
  
  # Compute acceptance probability (symmetric proposal)
  target_current <- dmvnorm(samples[t - 1, ], mean = mean_target, sigma = cov_target)
  target_candidate <- dmvnorm(candidate, mean = mean_target, sigma = cov_target)
  alpha <- min(1, target_candidate / target_current)
  
  # Accept or reject
  if (runif(1) < alpha) {
    samples[t, ] <- candidate
    accept_count <- accept_count + 1
  } else {
    samples[t, ] <- samples[t - 1, ]
  }
}

# ---- Acceptance Rate ----
accept_rate <- accept_count / (N - 1)
cat("Acceptance rate:", round(accept_rate, 4), "\n")

# ---- Plotting ----
# 1. Contour of Target Distribution
x1 <- seq(-3, 3, length.out = 100)
x2 <- seq(-3, 3, length.out = 100)
grid <- expand.grid(x1, x2)
z <- matrix(dmvnorm(grid, mean = mean_target, sigma = cov_target), 100, 100)

contour(x1, x2, z, nlevels = 10, col = "indianred1", lwd = 2,
        xlab = expression(theta[1]), ylab = expression(theta[2]),
        main = "RW-MH Samples with Target Contours")

# 2. Overlay samples (last 1000 for clarity)
points(samples[(N-999):N, 1], samples[(N-999):N, 2], 
       col = rgb(0, 0, 1, 0.3), pch = 19, cex = 0.5)

legend("topright", legend = c("Target Density", "RW-MH Samples"),
       col = c("indianred1", rgb(0, 0, 1, 0.3)), lwd = c(2,NA), pch = c(NA, 19), bty = "n")
