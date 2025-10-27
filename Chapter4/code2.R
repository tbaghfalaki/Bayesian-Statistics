# Load required library
library(ggplot2)

# -------------------------------
# Location prior (uniform)
theta_loc <- seq(-5, 5, length.out = 100)
pi_loc <- rep(1, length(theta_loc))

# Scale prior (1/theta)
theta_scale <- seq(0.01, 5, length.out = 100)
pi_scale <- 1/theta_scale

# Combine into one data frame for ggplot2
prior_df <- data.frame(
  theta = c(theta_loc, theta_scale),
  prior = c(pi_loc, pi_scale),
  type = factor(c(rep("Location", length(theta_loc)),
                  rep("Scale", length(theta_scale))))
)

# Plot both priors in one panel using facets
ggplot(prior_df, aes(x=theta, y=prior, color=type)) +
  geom_line(size=1.2) +
  facet_wrap(~type, scales = "free") +
  labs(title="Noninformative Priors for Location and Scale Parameters",
       x=expression(theta), y="Prior density") +
  theme_minimal() +
  theme(legend.position="none")




# Load required library
library(ggplot2)

# -------------------------------
# Location prior (uniform)
theta_loc <- seq(-5, 5, length.out = 100)
pi_loc <- rep(1, length(theta_loc))

# Scale prior (1/theta)
theta_scale <- seq(0.01, 5, length.out = 100)
pi_scale <- 1/theta_scale

# Combine into one data frame for ggplot2
prior_df <- data.frame(
  theta = c(theta_loc, theta_scale),
  prior = c(pi_loc, pi_scale),
  type = factor(c(rep("Location", length(theta_loc)),
                  rep("Scale", length(theta_scale))))
)

# Plot both priors in one panel using facets, without axis numbers
ggplot(prior_df, aes(x=theta, y=prior, color=type)) +
  geom_line(size=1.2) +
  facet_wrap(~type, scales = "free") +
  labs(title="Noninformative Priors for Location and Scale Parameters",
       x=expression(theta), y="Prior density") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),    # Remove axis numbers
    axis.ticks = element_blank()    # Remove axis ticks as well
  )

#######################################################################





