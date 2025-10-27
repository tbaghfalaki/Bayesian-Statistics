set.seed(123)

# Observed data
x <- c(5,6,7,5,6,6,7,5,6,6)
n <- length(x)
x_bar <- mean(x)
sigma2 <- 1
sigma <- sqrt(sigma2)

# Hyperpriors
mu0 <- 0
sigma_mu <- 100
scale_tau <- 5  # Half-Cauchy scale

# Number of Monte Carlo samples
n_sims <- 10000

# 1. Sample tau from Half-Cauchy
tau <- abs(rcauchy(n_sims, location = 0, scale = scale_tau))
tau2 <- tau^2

# 2. Sample mu from posterior conditional on tau and the data (marginalising theta)
# Xbar | mu, tau2 ~ N(mu, tau2 + sigma2/n)
mu_post_var <- 1 / ( 1/(tau2 + sigma2/n) + 1/sigma_mu^2 )
mu_post_mean <- mu_post_var * ( x_bar/(tau2 + sigma2/n) + mu0/sigma_mu^2 )
mu <- rnorm(n_sims, mean = mu_post_mean, sd = sqrt(mu_post_var))

# 3. Sample theta from its posterior given mu and tau2
theta_post_var <- 1 / ( n/sigma2 + 1/tau2 )
theta_post_mean <- theta_post_var * ( n*x_bar/sigma2 + mu/tau2 )
theta_samples <- rnorm(n_sims, mean = theta_post_mean, sd = sqrt(theta_post_var))

# 4. Posterior predictive samples
x_tilde_samples <- rnorm(n_sims, mean = theta_samples, sd = sigma)

# Posterior summaries
cat("Posterior mean of theta:         ", mean(theta_samples), "\n")
cat("Posterior SD of theta:           ", sd(theta_samples), "\n")
cat("Posterior predictive mean of ~X:  ", mean(x_tilde_samples), "\n")
cat("Posterior predictive SD of ~X:    ", sd(x_tilde_samples), "\n")
cat("Posterior mean of tau:           ", mean(tau), "\n")
cat("Posterior SD of tau:             ", sd(tau), "\n")

# Plot: single histogram + vertical lines + legend
hist(x_tilde_samples, breaks = 50,
     main = "Posterior Predictive Distribution",
     xlab = expression(tilde(X)),
     col = "skyblue", border = "white",
     xlim = range(c(x, x_tilde_samples)))

# add lines:
#abline(v = mean(x_tilde_samples), col = "red", lwd = 2)            # posterior predictive mean
#abline(v = mean(theta_samples), col = "darkgreen", lwd = 2, lty = 2) # posterior mean of theta
abline(v = x_bar, col = "red", lwd = 1.5)               # observed sample mean


