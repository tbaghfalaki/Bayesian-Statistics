# Load required packages
library(ggplot2)
library(extraDistr)

# Sequence of tau values
tau <- seq(0.01, 20, length.out = 1000)  # start at 0.01 to avoid division by zero

# 1. Half-Cauchy prior (scale = 5)
half_cauchy <- 2 / (pi * 5 * (1 + (tau/5)^2))

# 2. Inverse-Gamma prior on tau^2 (a = 0.1, b = 0.1) and transform to tau
inv_gamma <- dinvgamma(tau , a = 0.01, b = 0.01)    # Jacobian for tau^2 -> tau

# 3. Uniform prior (0 to 20)
uniform <- dunif(tau, min = 0, max = 20)

# 4. Non-informative prior proportional to 1/tau
noninformative <- 1 / tau
noninformative <- noninformative / sum(noninformative)  # normalize for plotting

# Combine into a data frame for plotting
prior_df <- data.frame(
  tau = rep(tau, 4),
  density = c(half_cauchy, inv_gamma, uniform, noninformative),
  prior = rep(c("Half-Cauchy", "Inverse-Gamma", "Uniform", "1/tau (non-informative)"),
              each = length(tau))
)

# Plot priors
ggplot(prior_df, aes(x = tau, y = density, color = prior)) +
  geom_line(size = 1.2) +
  labs(title = "Comparison of Priors for Tau",
       x = expression(tau),
       y = "Density") +
  theme_minimal()



# Plot priors with expression title
ggplot(prior_df, aes(x = tau, y = density, color = prior)) +
  geom_line(size = 1.2) +
  labs(title = expression("Comparison of Priors for " * tau),
       x = expression(tau),
       y = "Density") +
  theme_minimal()

