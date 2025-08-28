# Set the rate parameter
lambda <- 2

# --- Generate one sample using inverse transform sampling ---
u_single <- runif(1)
x_single <- -log(1 - u_single) / lambda

cat("Single uniform sample u:", u_single, "\n")
cat("Single exponential sample x:", x_single, "\n\n")

# --- Generate multiple samples ---
n <- 1000
u <- runif(n)
x <- -log(1 - u) / lambda  # Inverse CDF transformation

# --- Plot histogram with theoretical density ---
hist(x, breaks = 40, probability = TRUE, col = "lightblue",
     main = "Inverse Transform Sampling from Exp(2)",
     xlab = "x")

# Overlay the theoretical density function
curve(dexp(x, rate = lambda), col = "red", lwd = 2, add = TRUE)

# Add legend
legend("topright", legend = c("Sampled Histogram", "Theoretical Density"),
       col = c("lightblue", "red"), lwd = 2, bty = "n")
