library(plotly)

# Define grid for mu and sigma
mu <- seq(-5, 5, length.out = 100)
sigma <- seq(0.01, 5, length.out = 100)

# Create grid
grid <- expand.grid(mu = mu, sigma = sigma)

# Joint prior: π(μ, σ) = 1/σ
grid$prior <- 1 / grid$sigma

# Convert to matrix for plotly
prior_matrix <- matrix(grid$prior, nrow = length(mu), ncol = length(sigma))

# 3D plot using plotly with Greek symbols and no tick labels
fig <- plot_ly(x = ~mu, y = ~sigma, z = ~prior_matrix) %>%
  add_surface() %>%
  layout(
    title = "Noninformative Prior: Location–Scale",
    scene = list(
      xaxis = list(title = list(text = "\u03B8"), showticklabels = FALSE), # θ
      yaxis = list(title = list(text = "\u03C3"), showticklabels = FALSE), # σ
      zaxis = list(title = list(text = "\u03C0(\u03B8, \u03C3)"), showticklabels = FALSE) # π(θ, σ)
    )
  )

fig
