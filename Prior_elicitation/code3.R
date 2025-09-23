library(plotly)

# Define grid for mu and sigma
mu <- seq(-5, 5, length.out = 100)
sigma <- seq(0.01, 5, length.out = 100)

# Create grid
grid <- expand.grid(mu = mu, sigma = sigma)

# Joint prior: pi(mu, sigma) = 1/sigma
grid$prior <- 1 / grid$sigma

# Convert to matrix for plotly
prior_matrix <- matrix(grid$prior, nrow = length(mu), ncol = length(sigma))

# 3D plot using plotly with expressions and no numbers
fig <- plot_ly(x = ~mu, y = ~sigma, z = ~prior_matrix) %>%
  add_surface() %>%
  layout(
    title = "Noninformative Prior: Location-Scale",
    scene = list(
      xaxis = list(title = list(text = "μ"), showticklabels = FALSE),
      yaxis = list(title = list(text = "σ"), showticklabels = FALSE),
      zaxis = list(title = list(text = "π(μ, σ)"), showticklabels = FALSE)
    )
  )

fig
