
# Load the LearnBayes package
library(LearnBayes)

# Define midpoints for the prior
midpt <- seq(0, 1, by = 0.1)

# Assign subjective weights to each midpoint
prior <- c(0, 1, 2, 3, 4, 1, 2, 1, 0, 0)

# Normalize the prior weights to sum to 1
prior <- prior / sum(prior)

# Plot the prior density using the histogram method with color
curve(histprior(x, midpt, prior), 
      from = 0, to = 1,
      xlab = expression(paste(theta)),
      ylab = "Prior density",
      ylim = c(0, 0.3),
      main = "Prior Density via Histogram Method",
      col = "steelblue",      # Color of the curve
      lwd = 2)                # Line width
