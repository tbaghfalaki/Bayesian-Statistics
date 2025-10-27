# Load required library
library(mvtnorm)

# Target distribution parameters
mean_target <- c(0, 0)
rho <- -0.6
cov_target <- matrix(c(1, rho, rho, 1), nrow = 2)

# Gibbs sampling function
gibbs_sampler <- function(n_samples = 5000, rho = -0.6) {
  samples <- matrix(0, nrow = n_samples, ncol = 2)
  theta <- c(0, 0)  # initial value
  sd_cond <- sqrt(1 - rho^2)
  
  for (n in 2:n_samples) {
    # Sample theta_1 | theta_2
    theta[1] <- rnorm(1, mean = rho * theta[2], sd = sd_cond)
    
    # Sample theta_2 | theta_1
    theta[2] <- rnorm(1, mean = rho * theta[1], sd = sd_cond)
    
    samples[n, ] <- theta
  }
  
  return(samples)
}

# Run Gibbs sampler
set.seed(123)
gibbs_samples <- gibbs_sampler()

# Contour plotting function
plot_contour <- function(mean_vec, cov_matrix, samples = NULL) {
  x1 <- seq(-3, 3, 0.1)
  x2 <- seq(-3, 3, 0.1)
  grid <- expand.grid(x1, x2)
  z <- matrix(dmvnorm(grid, mean = mean_vec, sigma = cov_matrix), nrow = length(x1))
  
  contour(x1, x2, z, nlevels = 10, col = "indianred1", lwd = 2,
          xlab = expression(theta[1]), ylab = expression(theta[2]),
          main = "Target Contours and Gibbs Samples")
  
  if (!is.null(samples)) {
    points(samples, pch = 19, col = rgb(0, 0, 1, 0.3), cex = 0.5)
  }
}

# Plot contours with last 1000 samples
plot_contour(mean_target, cov_target, gibbs_samples[4001:5000, ])
