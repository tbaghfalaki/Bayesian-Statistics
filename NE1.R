# Load necessary library
library(ggplot2)

# Given parameters
alpha_prior <- 1
beta_prior <- 1
n <- 10
sum_x <- 6

# Compute posterior parameters
alpha_post <- alpha_prior + sum_x  # α* = α + sum(x)
beta_post <- beta_prior + n        # β* = β + n

# Define credibility level
credibility_level <- 0.95
alpha_half <- (1 - credibility_level)

# Define prior density function
prior_density <- function(theta) {
  (beta_prior^alpha_prior / gamma(alpha_prior)) * theta^(alpha_prior - 1) * exp(-beta_prior * theta)
}

# Define posterior density function
posterior_density <- function(theta) {
  (beta_post^alpha_post / gamma(alpha_post)) * theta^(alpha_post - 1) * exp(-beta_post * theta)
}

# Find k: numerical search for the threshold
k_values <- seq(0, posterior_density(alpha_post / beta_post), length.out = 1000)  # Possible values of k
intervals <- sapply(k_values, function(k) {
  theta_range <- seq(0, 2, length.out = 10000)
  density_values <- posterior_density(theta_range)
  selected <- theta_range[density_values > k]
  if (length(selected) > 1) {
    return(c(min(selected), max(selected)))
  } else {
    return(c(NA, NA))
  }
})

# Find the shortest interval containing 95% probability
valid_intervals <- intervals[, !is.na(intervals[1, ])]
probs <- apply(valid_intervals, 2, function(bounds) {
  integrate(posterior_density, lower = bounds[1], upper = bounds[2])$value
})

# Choose the k corresponding to the shortest valid HPD interval
best_index <- which.min(abs(probs - credibility_level))
hpd_L <- valid_intervals[1, best_index]
hpd_U <- valid_intervals[2, best_index]
k_best <- k_values[best_index]

# Compute equal-tailed credible interval
theta_L <- qgamma(alpha_half / 2, shape = alpha_post, rate = beta_post)
theta_U <- qgamma(1 - alpha_half / 2, shape = alpha_post, rate = beta_post)

# Print results
cat("Equal-tailed credible interval: (", theta_L, ",", theta_U, ")\n")
cat("HPD interval: (", hpd_L, ",", hpd_U, ")\n")
cat("Threshold k:", k_best, "\n")

# Generate density plots
theta_vals <- seq(0, 2, length.out = 1000)
posterior_vals <- posterior_density(theta_vals)
prior_vals <- prior_density(theta_vals)

df_posterior <- data.frame(theta_vals, posterior_vals)
df_prior <- data.frame(theta_vals, prior_vals)

# Plot 1: Prior and Posterior
ggplot() +
  geom_line(data = df_prior, aes(x = theta_vals, y = prior_vals), color = "black", linetype = "dashed", size = 1.2) +
  geom_line(data = df_posterior, aes(x = theta_vals, y = posterior_vals), color = "darkorchid", size = 1.2) +
  labs(title = "Prior and Posterior Distributions",
       x = "θ", 
       y = "Density") +
  theme_minimal() +
  annotate("text", x = 1.2, y = max(prior_vals) * 0.4, label = "Prior", color = "black", hjust = -0.1) +
  annotate("text", x = 0.5, y = max(posterior_vals) * 0.9, label = "Posterior", color = "darkorchid", hjust = -0.1)

# Plot 2: Posterior with HPD and Equal-Tailed Intervals
ggplot(df_posterior, aes(x = theta_vals, y = posterior_vals)) +
  geom_line(color = "darkorchid", size = 1.2) +
  geom_vline(xintercept = c(theta_L, theta_U), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = c(hpd_L, hpd_U), linetype = "dashed", color = "green", size = 1) +
  labs(title = "Posterior Distribution with HPD and Equal-Tailed Intervals",
       x = "θ", 
       y = "Density") +
  theme_minimal() +
  annotate("text", x = theta_L, y = max(posterior_vals) * 0.8, label = "Lower ET", color = "red", hjust = -0.1) +
  annotate("text", x = theta_U, y = max(posterior_vals) * 0.8, label = "Upper ET", color = "red", hjust = 1.1) +
  annotate("text", x = hpd_L, y = max(posterior_vals) * 0.6, label = "Lower HPD", color = "green", hjust = -0.1) +
  annotate("text", x = hpd_U, y = max(posterior_vals) * 0.6, label = "Upper HPD", color = "green", hjust = 1.1)
