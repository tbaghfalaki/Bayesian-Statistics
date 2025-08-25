# Load required library
library(ggplot2)
library(gridExtra)

# Define parameters
n <- 10
x_bar <- 5
sigma2 <- 1
sigma <- sqrt(sigma2)
mu0 <- 0
sigma02 <- 100
sigma0 <- sqrt(sigma02)

# Compute posterior parameters
posterior_mean <- (n * x_bar / sigma2 + mu0 / sigma02) / (n / sigma2 + 1 / sigma02)
posterior_var <- 1 / (n / sigma2 + 1 / sigma02)
posterior_sd <- sqrt(posterior_var)

# Compute posterior predictive distribution parameters
predictive_mean <- posterior_mean
predictive_var <- posterior_var + sigma2
predictive_sd <- sqrt(predictive_var)

# Define x values for plotting
theta_vals <- seq(1, 9, length = 1000)

# Compute densities
prior_density <- dnorm(theta_vals, mean = mu0, sd = sigma0)
posterior_density <- dnorm(theta_vals, mean = posterior_mean, sd = posterior_sd)


# Create data frames for plotting
prior_posterior_data <- data.frame(
  x = rep(theta_vals, 2),
  density = c(prior_density, posterior_density),
  Distribution = rep(c("Prior", "Posterior"), each = length(theta_vals))
)


# Plot prior and posterior distributions
ggplot(prior_posterior_data, aes(x = x, y = density, color = Distribution, linetype = Distribution)) +
  geom_line(size = 1) +
  labs(title = "Prior and Posterior Distributions",
       x = expression(theta),
       y = "Density") +
  theme_minimal() +
  scale_color_manual(values = c("Prior" = "#1f77b4", "Posterior" = "#d62728")) +
  scale_linetype_manual(values = c("Prior" = "solid", "Posterior" = "solid"))

x_vals <- seq(1, 9, length = 1000)


# Compute Bayesian posterior predictive distribution
predictive_density <- dnorm(x_vals, mean = predictive_mean, sd = predictive_sd)

predictive_data <- data.frame(
  x = x_vals,
  density = predictive_density,
  Distribution = "Posterior Predictive"
)


# Plot posterior predictive distribution
ggplot(predictive_data, aes(x = x, y = density, color = Distribution)) +
  geom_line(size = 1) +
  labs(title = "Posterior Predictive Distribution",
       x = expression(x[n+1]),
       y = "Density") +
  theme_minimal() +
  scale_color_manual(values = c("Posterior Predictive" = "#2ca02c"))