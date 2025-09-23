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

# 2. Sample mu from its conditional posterior given tau
mu_post_sd <- sqrt(1 / (1/tau2 + 1/sigma_mu^2))
mu_post_mean <- mu_post_sd^2 * (0/sigma_mu^2 + x_bar/tau2)
mu <- rnorm(n_sims, mean = mu_post_mean, sd = mu_post_sd)

# 3. Sample theta from its conditional posterior given mu and tau2
theta_post_sd <- sqrt(1 / (n/sigma2 + 1/tau2))
theta_post_mean <- theta_post_sd^2 * (n*x_bar/sigma2 + mu/tau2)
theta_samples <- rnorm(n_sims, mean = theta_post_mean, sd = theta_post_sd)

# 4. Posterior predictive samples
x_tilde_samples <- rnorm(n_sims, mean = theta_samples, sd = sigma)

# Posterior summaries
cat("Posterior mean of theta:", mean(theta_samples), "\n")
cat("Posterior SD of theta:", sd(theta_samples), "\n")
cat("Posterior predictive mean of ~X:", mean(x_tilde_samples), "\n")
cat("Posterior predictive SD of ~X:", sd(x_tilde_samples), "\n")

# Histogram of posterior predictive samples
hist(x_tilde_samples, breaks = 50,
     main = "Posterior Predictive Distribution",
     xlab = expression(tilde(X)),
     col = "skyblue", border = "white",
     xlim = range(c(x, x_tilde_samples)))

# Add posterior predictive mean
abline(v = mean(x_tilde_samples), col = "red", lwd = 2)
