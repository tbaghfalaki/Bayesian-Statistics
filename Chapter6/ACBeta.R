rm(list=ls())
ACBeta <- function(n) {
  u <- runif(n)
  x <- runif(n)
  f <- 6 * x * (1 - x)
  accepted_samples <- x[u < 2 * f / 3]
  return(accepted_samples)
}

Y <- ACBeta(10000)

hist(Y, breaks = 30, probability = TRUE,
     main = "Histogram and Estimated Density",
     xlab = "x", col = "lightblue", border = "white")

# Overlay true Beta(2,2) density
curve(dbeta(x, 2, 2), add = TRUE, col = "red", lwd = 2)


hist(Y, breaks = 30, probability = TRUE,
     main = "Histogram and True Density",
     xlab = "x", col = "lightblue", border = "white")

# Overlay true Beta(2,2) density
curve(dbeta(x, 2, 2), add = TRUE, col = "red", lwd = 2)

# Add legend
legend("topright",
       legend = c("Histogram of Samples", "True Beta(2,2) Density"),
       fill = c("lightblue", "red"),  # 'fill' for histogram
       border = c("white", NA),
       #lty = c(NA, 1),              # line type for density
       col = c("lightblue", "red"),
       #lwd = c(2, 2),
       bty = "n")
