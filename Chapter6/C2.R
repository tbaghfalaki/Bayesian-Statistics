set.seed(123)

# --- Observed Data ---
n <- 100
theta_true <- 0.7
y <- rbinom(n, 1, theta_true)

# --- Prior Parameters ---
alpha <- 2
beta <- 2

# --- Derived Posterior Parameters ---
alpha_star <- alpha - 1 + sum(y)          # α*
beta_star  <- beta - 1 + n - sum(y)       # β*

# --- HMC Parameters ---
epsilon <- 0.05   # step size
L <- 20           # number of leapfrog steps
S <- 5000         # number of HMC samples

# --- Logistic transform ---
theta_from_z <- function(z) 1 / (1 + exp(-z))
dtheta_dz <- function(z) theta_from_z(z) * (1 - theta_from_z(z))  # derivative

# --- Potential energy in z-space ---
U <- function(z) {
  theta <- theta_from_z(z)
  - alpha_star * log(theta) -
    beta_star * log(1 - theta) -
    log(theta * (1 - theta))      # log-Jacobian term
}

# --- Gradient of potential energy wrt z ---
grad_U <- function(z) {
  theta <- theta_from_z(z)
  dtheta <- dtheta_dz(z)
  
  # chain rule
  (-alpha_star / theta + beta_star / (1 - theta)) * dtheta - 
    (1 - 2 * theta)                # derivative of log(theta(1-theta))
}

# --- HMC Sampling ---
samples <- numeric(S)
z <- qlogis(0.5)   # initial z corresponds to theta = 0.5

for (s in 1:S) {
  
  current_z <- z
  p <- rnorm(1)           # sample momentum
  current_p <- p          # save for Hamiltonian
  
  # Leapfrog integration
  p <- p - (epsilon / 2) * grad_U(z)
  
  for (l in 1:L) {
    z <- z + epsilon * p
    if (l != L) p <- p - epsilon * grad_U(z)
  }
  
  p <- p - (epsilon / 2) * grad_U(z)
  p <- -p    # momentum flip
  
  # Hamiltonians
  current_H <- U(current_z) + (current_p^2) / 2
  proposed_H <- U(z) + (p^2) / 2
  
  # Accept/reject
  accept_prob <- exp(current_H - proposed_H)
  
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
     main = "Full Trace",
     ylab = expression(theta), xlab = "Iteration")

dens <- density(samples)
plot(dens, col="darkgreen",
     main="Posterior Density",
     xlab=expression(theta), ylab="Density")
abline(v=theta_true, lty=2, col="red", lwd=2)
