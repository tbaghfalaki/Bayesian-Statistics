library(mvtnorm)

# Function to plot contour lines of two bivariate normal densities
plot_bivariate_normals <- function(mean1, cov1, mean2, cov2) {
  x_seq <- seq(-4, 4, 0.1)
  y_seq <- seq(-4, 4, 0.1)
  grid <- expand.grid(x_seq, y_seq)
  colnames(grid) <- c("x", "y")
  
  dens_target <- matrix(dmvnorm(grid, mean = mean1, sigma = cov1), length(x_seq), length(y_seq))
  dens_proposal <- matrix(dmvnorm(grid, mean = mean2, sigma = cov2), length(x_seq), length(y_seq))
  
  contour(x_seq, y_seq, dens_target, nlevels = 10, col = "indianred1", lwd = 2,
          xlab = expression(theta[1]), ylab = expression(theta[2]),
          main = "Target and Proposal Contours")
  contour(x_seq, y_seq, dens_proposal, nlevels = 10, col = "lightblue", lwd = 2, add = TRUE)
  
  legend("topleft", legend = c("Target density", "Proposal density"),
         col = c("indianred1", "lightblue"), lwd = 2, bty = "n")
}

# Function to plot target contour with samples overlay
plot_target_with_samples <- function(mean_target, cov_target, samples, nlevels = 10) {
  x_seq <- seq(-4, 4, 0.1)
  y_seq <- seq(-4, 4, 0.1)
  grid <- expand.grid(x_seq, y_seq)
  colnames(grid) <- c("x", "y")
  
  dens_target <- matrix(dmvnorm(grid, mean = mean_target, sigma = cov_target), length(x_seq), length(y_seq))
  
  contour(x_seq, y_seq, dens_target, nlevels = nlevels, col = "indianred1", lwd = 2,
          xlab = expression(theta[1]), ylab = expression(theta[2]),
          main = "Target Density with MH Samples")
  
  points(samples[,1], samples[,2], pch = 19, col = rgb(0, 0, 1, 0.3), cex = 0.5)
  
  legend("topleft", legend = c("Target density", "MH samples"),
         col = c("indianred1", rgb(0, 0, 1, 0.3)), pch = c(NA, 19), lwd = c(2, NA), bty = "n")
}

# Parameters with correlation -0.6 in target covariance
mean_target <- c(0, 0)
cov_target <- matrix(c(1, -0.6, -0.6, 1), 2, 2)

# Proposal distribution (independent standard normal)
mean_proposal <- c(0, 0)
cov_proposal <- diag(2)

# MH sampling setup
n_samples <- 2000
samples <- matrix(0, n_samples, 2)
samples[1, ] <- c(0, 0)
accept_count <- 0

for (t in 2:n_samples) {
  candidate <- rmvnorm(1, mean_proposal, cov_proposal)
  target_candidate <- dmvnorm(candidate, mean_target, cov_target)
  target_current <- dmvnorm(samples[t-1, ], mean_target, cov_target)
  alpha <- min(1, target_candidate / target_current)
  
  if (runif(1) < alpha) {
    samples[t, ] <- candidate
    accept_count <- accept_count + 1
  } else {
    samples[t, ] <- samples[t-1, ]
  }
}

accept_rate <- accept_count / (n_samples - 1)
cat("Acceptance rate:", round(accept_rate, 4), "\n")

# Plotting 2 panels
par(mfrow = c(1, 2), mar = c(4, 4, 4, 2))

# Panel 1: Target and proposal contours
plot_bivariate_normals(mean_target, cov_target, mean_proposal, cov_proposal)

# Panel 2: Target contour + last 1000 MH samples
plot_target_with_samples(mean_target, cov_target, samples[(n_samples-999):n_samples, ])
