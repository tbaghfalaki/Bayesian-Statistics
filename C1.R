set.seed(123)

# --- Parameters ---
n <- 100
theta_true <- 1
sigma2 <- 1
tau2 <- 100
epsilon <- 0.1   # Step size for leapfrog
L <- 20          # Number of leapfrog steps
S <- 2000        # Number of HMC samples

# --- Simulate data ---
y <- rnorm(n, mean = theta_true, sd = sqrt(sigma2))

# --- Potential energy function ---
U <- function(theta) {
  (theta^2) / (2 * tau2) + sum((y - theta)^2) / (2 * sigma2)
}

# --- Gradient of potential energy ---
grad_U <- function(theta) {
  theta / tau2 + (n * theta - sum(y)) / sigma2
}

# --- HMC sampling ---
samples <- numeric(S)
theta <- 2  # Initial parameter value

for (s in 1:S) {
  current_theta <- theta
  p <- rnorm(1)        # Sample momentum
  current_p <- p
  
  # Leapfrog integration
  p <- p - (epsilon / 2) * grad_U(theta)
  for (i in 1:L) {
    theta <- theta + epsilon * p
    if (i != L) p <- p - epsilon * grad_U(theta)
  }
  p <- p - (epsilon / 2) * grad_U(theta)
  p <- -p  # Momentum flip for reversibility
  
  # Calculate Hamiltonians
  current_H <- U(current_theta) + (current_p^2) / 2
  proposed_H <- U(theta) + (p^2) / 2
  
  # Accept/reject step
  if (runif(1) < exp(current_H - proposed_H)) {
    samples[s] <- theta
  } else {
    theta <- current_theta
    samples[s] <- theta
  }
}

# --- Posterior summaries ---
posterior_mean <- mean(samples)
posterior_sd <- sd(samples)

cat(sprintf("Posterior mean estimate: %.4f\n", posterior_mean))
cat(sprintf("Posterior standard deviation: %.4f\n", posterior_sd))

# --- Plots setup ---
first_jumps <- samples[1:50]

par(mfrow = c(3, 1), mar = c(4, 4, 3, 2))  # 3 rows, 1 col with margins

# 1) Zoomed trace plot (first 50 samples)
plot(first_jumps, type = "o", col = "dodgerblue3", pch = 16,
     main = "HMC Initial Jumps (First 50 Iterations)",
     xlab = "Iteration", ylab = expression(theta),
     ylim = range(c(first_jumps, theta_true)))
abline(h = theta_true, col = "red", lty = 2, lwd = 2)
legend("bottomright", legend = c("Sampled θ", "True θ"),
       col = c("dodgerblue3", "red"), lty = c(1, 2), pch = c(16, NA),
       bg = "white")

# 2) Full trace plot
plot(samples, type = "l", col = "dodgerblue3",
     main = "HMC Trace Plot for θ",
     xlab = "Iteration", ylab = expression(theta))
abline(h = theta_true, col = "red", lty = 2, lwd = 2)
legend("bottomright", legend = c("Sampled θ", "True θ"),
       col = c("dodgerblue3", "red"), lty = c(1, 2), bg = "white")

# 3) Posterior density plot
dens <- density(samples)
plot(dens, main = "Posterior Density of θ from HMC",
     xlab = expression(theta), ylab = "Density",
     col = "forestgreen", lwd = 2)
abline(v = theta_true, col = "red", lty = 2, lwd = 2)
legend("topright", legend = c("Estimated density", "True θ"),
       col = c("forestgreen", "red"), lty = c(1, 2), lwd = 2, bg = "white")

# Add posterior summaries text inside density plot
text(x = min(samples), y = max(dens$y) * 0.7,
     labels = sprintf("Posterior mean = %.3f\nPosterior SD = %.3f",
                      posterior_mean, posterior_sd),
     adj = c(0, 0), col = "black", cex = 0.9)
