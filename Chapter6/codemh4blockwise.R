library(mvtnorm)

# Target distribution parameters (bivariate normal with correlation)
mean_target <- c(0, 0)
cov_target <- matrix(c(1, -0.6, -0.6, 1), nrow = 2)
inv_cov <- solve(cov_target)

# Unnormalized log target density function
log_target_density <- function(theta) {
  -0.5 * t(theta) %*% inv_cov %*% theta
}

# Blockwise Metropolis-Hastings sampler
blockwise_mh <- function(n_samples = 5000, sigma1 = 0.4, sigma2 = 0.4) {
  samples <- matrix(0, nrow = n_samples, ncol = 2)
  accept_counts <- c(0, 0)
  
  # Initialize chain
  theta <- c(0, 0)
  samples[1, ] <- theta
  
  for (n in 2:n_samples) {
    for (j in 1:2) {
      # Propose new value for component j
      theta_prop <- theta
      proposal_sd <- ifelse(j == 1, sigma1, sigma2)
      theta_prop[j] <- rnorm(1, mean = theta[j], sd = proposal_sd)
      
      # Compute acceptance probability (log scale)
      log_alpha <- log_target_density(theta_prop) - log_target_density(theta)
      alpha <- min(1, exp(log_alpha))
      
      # Accept or reject
      if (runif(1) < alpha) {
        theta <- theta_prop
        accept_counts[j] <- accept_counts[j] + 1
      }
    }
    samples[n, ] <- theta
  }
  
  list(samples = samples, acceptance_rate = accept_counts / (n_samples - 1))
}

# Run the sampler
set.seed(123)
result <- blockwise_mh(n_samples = 5000, sigma1 = 0.4, sigma2 = 0.4)

# Acceptance rates
cat("Acceptance rates:\n")
cat(sprintf("theta_1: %.3f\n", result$acceptance_rate[1]))
cat(sprintf("theta_2: %.3f\n", result$acceptance_rate[2]))

# Plotting function for contours and samples
plot_contours_and_samples <- function(mean_vec, cov_mat, samples, burnin = 1000) {
  x <- seq(-3, 3, length.out = 100)
  y <- seq(-3, 3, length.out = 100)
  grid <- expand.grid(x, y)
  z <- matrix(dmvnorm(grid, mean = mean_vec, sigma = cov_mat), nrow = length(x))
  
  contour(x, y, z, nlevels = 10, col = "red", lwd = 2,
          xlab = expression(theta[1]), ylab = expression(theta[2]),
          main = "Target Distribution Contours and Blockwise MH Samples")
  
  points(samples[(burnin + 1):nrow(samples), ], pch = 19, col = rgb(0, 0, 1, 0.3), cex = 0.5)
}

# Plot the results (last 4000 samples after burn-in)
plot_contours_and_samples(mean_target, cov_target, result$samples, burnin = 4000)
