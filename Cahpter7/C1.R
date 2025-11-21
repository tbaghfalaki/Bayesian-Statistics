### ---- Variational Bayes for Normal Model (Example 2) ----
rm(list=ls())
set.seed(9)

# --- Generate data ---
n <- 100
theta_true <- 2
sigma2_true <- 1
y <- rnorm(n, mean = theta_true, sd = sqrt(sigma2_true))

# --- Prior parameters ---
tau2 <- 100         # prior variance for theta
a0 <- 0.01            # IG prior shape
b0 <- 0.01            # IG prior scale

# --- Initialize VB parameters ---
mu <- 0
lambda <- 1
tilde_a <- a0
tilde_b <- b0

max_iter <- 2000
tol <- 1e-8

# --- Coordinate Ascent Variational Inference ---
for (iter in 1:max_iter) {
  
  mu_old <- mu
  lambda_old <- lambda
  tb_old <- tilde_b
  
  # update q(theta)
  lambda <- n * (tilde_a / tilde_b) + 1/tau2
  mu <- ( (tilde_a / tilde_b) * sum(y) ) / lambda
  
  # update q(sigma^2)
  tilde_a <- a0 + n/2
  tilde_b <- b0 + 0.5 * sum( (y - mu)^2 + 1/lambda )
  
  # convergence check
  if (max(abs(mu - mu_old), abs(lambda - lambda_old), abs(tilde_b - tb_old)) < tol)
    break
}
print(iter) # convergence check stop here

# Posterior mean and SD analytically
theta_mean_vb <- mu
theta_sd_vb   <- sqrt(1/lambda)

sigma2_mean_vb <- tilde_b / (tilde_a - 1)
sigma2_sd_vb   <- sqrt(tilde_b^2 / ((tilde_a - 1)^2 * (tilde_a - 2)))

cat("VB posterior mean of theta:", theta_mean_vb, "\n")
cat("VB posterior SD of theta:", theta_sd_vb, "\n\n")

cat("VB posterior mean of sigma^2:", sigma2_mean_vb, "\n")
cat("VB posterior SD of sigma^2:", sigma2_sd_vb, "\n")





# --- Optional: Draw posterior samples from q* (if you want to plot later) ---
S <- 10000
theta_samples <- rnorm(S, mean = mu, sd = sqrt(1/lambda))
sigma2_samples <- 1 / rgamma(S, shape = tilde_a, rate = tilde_b)




library(ggplot2)

library(ggplot2)
library(dplyr)

# Prepare data in long format
df <- data.frame(
  value = c(theta_samples, sigma2_samples),
  parameter = rep(c("theta", "sigma2"), each = length(theta_samples))
)

# True values
true_vals <- data.frame(
  parameter = c("theta", "sigma2"),
  true = c(theta_true, sigma2_true)
)

# Plot with facets
p_facet <- ggplot(df, aes(x = value, fill = parameter)) +
  geom_density(alpha = 0.5) +
  geom_vline(data = true_vals, aes(xintercept = true, color = parameter),
             linetype = "dashed", size = 1) +
  facet_wrap(~parameter, scales = "free") +
  scale_fill_manual(values = c("skyblue", "yellow3")) +
  scale_color_manual(values = c("blue", "orange3")) +
  labs(title = "VB Posterior Densities of θ and σ²",
       x = "Value", y = "Density") +
  theme_bw() +
  theme(legend.position = "none")

print(p_facet)

