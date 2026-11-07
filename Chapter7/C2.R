set.seed(123)

# --- Observed Data ---
n <- 100
theta_true <- 0.7
y <- rbinom(n, 1, theta_true)

# --- Prior Parameters ---
alpha <- 2
beta <- 2

# --- Derived Posterior Parameters ---
alpha_star <- alpha - 1 + sum(y)
beta_star  <- beta - 1 + n - sum(y)

# --- HMC Parameters ---
epsilon <- 0.05
L <- 20
S <- 5000

# --- Logistic transform ---
theta_from_z <- function(z) 1 / (1 + exp(-z))
dtheta_dz <- function(z) theta_from_z(z) * (1 - theta_from_z(z))

# --- Potential energy in z-space ---
U <- function(z) {
  theta <- theta_from_z(z)
  -alpha_star * log(theta) - beta_star * log(1 - theta) - log(theta * (1 - theta))
}

# --- Gradient of potential energy wrt z ---
grad_U <- function(z) {
  theta <- theta_from_z(z)
  dtheta <- dtheta_dz(z)
  # Chain rule including Jacobian
  (-alpha_star / theta + beta_star / (1 - theta) - 1/theta + 1/(1-theta)) * dtheta
}

# --- HMC Sampling ---
samples <- numeric(S)
z <- qlogis(0.5)   # initial z

for (s in 1:S) {
  current_z <- z
  p <- rnorm(1)
  current_p <- p
  
  # Leapfrog
  p <- p - (epsilon/2) * grad_U(z)
  for (l in 1:L) {
    z <- z + epsilon * p
    if (l != L) p <- p - epsilon * grad_U(z)
  }
  p <- p - (epsilon/2) * grad_U(z)
  p <- -p   # momentum flip
  
  # Hamiltonian
  current_H <- U(current_z) + (current_p^2)/2
  proposed_H <- U(z) + (p^2)/2
  
  # Accept/reject
  accept_prob <- min(1, exp(current_H - proposed_H))
  if (runif(1) < accept_prob) {
    samples[s] <- theta_from_z(z)
  } else {
    z <- current_z
    samples[s] <- theta_from_z(z)
  }
}

# --- Posterior summaries ---
posterior_mean <- mean(samples)
posterior_sd <- sd(samples)

cat(sprintf("Posterior mean estimate: %.4f\n", posterior_mean))
cat(sprintf("Posterior standard deviation: %.4f\n", posterior_sd))

# --- Diagnostic plots ---
par(mfrow = c(2, 1))
plot(samples, type = "l", col = "blue",
     main = "Trace of θ",
     ylab = expression(theta), xlab = "Iteration")

dens <- density(samples)
plot(dens, col = "darkgreen",
     main = "Posterior Density of θ",
     xlab = expression(theta), ylab = "Density")
abline(v = theta_true, lty = 2, col = "red", lwd = 2)
