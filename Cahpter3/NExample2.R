# Load required library
library(ggplot2)
library(tibble)

# Define prior parameters
alpha <- 0.01
beta <- 0.01

# Observed data
X <- c(1.2, 0.8, 1.5, 1.0, 1.3)
n <- length(X)
S <- sum(X)

# Compute posterior parameters
alpha_post <- alpha + n
beta_post <- beta + S

# Define range for theta
theta_vals <- seq(0, 3, length = 1000)

# Compute densities
prior_density <- dgamma(theta_vals, shape = alpha, rate = beta)
posterior_density <- dgamma(theta_vals, shape = alpha_post, rate = beta_post)

# Create a data frame for plotting
prior_data <- data.frame(theta = theta_vals, density = prior_density, Distribution = "Prior")
posterior_data <- data.frame(theta = theta_vals, density = posterior_density, Distribution = "Posterior")

data <- rbind(prior_data, posterior_data)

# Plot prior and posterior distributions
ggplot(data, aes(x = theta, y = density, color = Distribution, linetype = Distribution)) +
  geom_line(size = 1) +
  labs(title = "Prior and Posterior Distributions",
       x = expression(theta),
       y = "Density") +
  theme_minimal() +
  scale_color_manual(values = c("Prior" = "blue", "Posterior" = "red")) +
  scale_linetype_manual(values = c("Prior" = "dashed", "Posterior" = "solid"))

# Posterior predictive distribution

# Generate sequence
x_new <- seq(0, 10, length.out = 10000)

# Compute posterior predictive samples
posterior_predictive <- (5.81^5.01) / (5.81 + x_new)^6.01

# Create a tibble (data frame) for ggplot
data <- tibble(x = x_new, y = posterior_predictive)

# Create the plot
ggplot(data, aes(x = x, y = y)) +
  geom_line(color = "blue",size = 1) +
  labs(title = "Posterior Predictive Distribution",
       x = expression(x[new]),
       y = "Density") +
  theme_minimal()