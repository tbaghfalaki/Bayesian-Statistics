# Load necessary libraries
library(ggplot2)
library(HDInterval)  # Package for HPD interval

# Given data
mu0 <- 0                      # Prior mean
sigma2 <- 4                   # Known variance of likelihood
n <- 10                       # Sample size
xbar <- 2.5                   # Sample mean
sigma0_sq_vals <- c(0.1, 1, 10, 100, 1000)  # Different prior variances
alpha <- 0.05                 # Significance level for 95% intervals

# Define theta values for plotting
theta_vals <- seq(-2, 8, length.out = 500)

# Function to compute posterior mean and variance
posterior_params <- function(sigma0_sq) {
  mu_n <- (mu0/sigma0_sq + n*xbar/sigma2) / (1/sigma0_sq + n/sigma2)
  sigma_n_sq <- 1 / (1/sigma0_sq + n/sigma2)
  return(data.frame(sigma0_sq, mu_n, sigma_n_sq))
}

# Compute posterior parameters for each sigma0_sq
results <- do.call(rbind, lapply(sigma0_sq_vals, posterior_params))

# Function to compute Equal-Tailed and HPD intervals using the HDInterval package
compute_intervals <- function(mu_n, sigma_n_sq) {
  sigma_n <- sqrt(sigma_n_sq)
  
  # Equal-Tailed Interval (ETI)
  eti_lower <- qnorm(alpha/2, mean = mu_n, sd = sigma_n)
  eti_upper <- qnorm(1 - alpha/2, mean = mu_n, sd = sigma_n)
  
  # HPD Interval using HDInterval package
  hpd_vals <- hdi(rnorm(100000, mean = mu_n, sd = sigma_n), credMass = 1 - alpha)
  
  return(data.frame(ETI_Lower = eti_lower, ETI_Upper = eti_upper,
                    HPD_Lower = hpd_vals[1], HPD_Upper = hpd_vals[2]))
}

# Compute intervals for each prior variance
credible_intervals <- do.call(rbind, 
                              mapply(compute_intervals, results$mu_n, results$sigma_n_sq, SIMPLIFY = FALSE))
credible_intervals$PriorVariance <- results$sigma0_sq

# Data for posterior distributions
posterior_densities <- data.frame()

for (i in 1:nrow(results)) {
  sigma0_sq <- results$sigma0_sq[i]
  mu_n <- results$mu_n[i]
  sigma_n <- sqrt(results$sigma_n_sq[i])
  
  # Compute posterior density
  posterior_dens <- dnorm(theta_vals, mean = mu_n, sd = sigma_n)
  
  # Store in dataframe
  posterior_densities <- rbind(posterior_densities,
                               data.frame(theta_vals, Density = posterior_dens, 
                                          PriorVariance = factor(sigma0_sq)))
}

# Merge posterior densities with credible intervals for faceted plot
plot_data <- merge(posterior_densities, credible_intervals, by = "PriorVariance")

# Plot Posterior Distributions in a Single Column with Fixed Scale
ggplot(plot_data, aes(x = theta_vals, y = Density)) +
  geom_line(size = 1, color = "blue") +
  geom_vline(aes(xintercept = ETI_Lower), linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = ETI_Upper), linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = HPD_Lower), linetype = "dotted", color = "black") +
  geom_vline(aes(xintercept = HPD_Upper), linetype = "dotted", color = "black") +
  facet_wrap(~PriorVariance, ncol = 1, scales = "fixed") +  # One column layout
  labs(title = "Posterior Distributions with 95% Credible Intervals",
       subtitle = "Dashed = Equal-Tailed, Dotted = HPD",
       x = "Theta (θ)", y = "Density") +
  theme_minimal()