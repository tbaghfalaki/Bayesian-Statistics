rm(list=ls())
# Improved comparison of Squared Error Loss and Stein Loss
par(mfrow = c(1, 2), mar = c(4, 5, 3, 5), oma = c(0, 0, 2, 0))

# Function to plot both losses for a given true value 'a'
plot_losses <- function(a, panel_label) {
  sigma2 <- seq(0.1, 4, length.out = 300)
  square_error_loss <- (sigma2 - a)^2
  stein_loss <- (a / sigma2) - 1 - log(a / sigma2)
  
  plot(sigma2, square_error_loss, type = "l", lwd = 2, col = "blue",
       ylim = c(0, max(c(square_error_loss, stein_loss))),
       xlab = expression(sigma^2), ylab = expression(L(sigma^2, a)),
       main = panel_label, cex.main = 1.2, cex.lab = 1.2)
  
  lines(sigma2, stein_loss, col = "red", lwd = 2)
  abline(v = a, lty = 2, col = "black")
}

# Plot for a = 1 and a = 2
plot_losses(a = 1, panel_label = "")
plot_losses(a = 2, panel_label = "")

# Shared legend
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")


legend("topright",bty = "n",
       legend = c(expression((sigma^2 - a)^2),
                  expression( a/sigma^2 - 1 - log(a/sigma^2)),
                  expression(True~value~of~a)),
       col = c("blue", "red", "black"), lwd = c(2, 2, 1), lty = c(1, 1, 2))
