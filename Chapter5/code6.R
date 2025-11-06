library(R2jags)

# Observed data
x <- c(5, 6, 7, 5, 6, 6, 7, 5, 6, 6)
n <- length(x)
sigma2 <- 1  # known
sigma <- sqrt(sigma2)

# JAGS model with Inverse-Gamma prior for tau^2
model_string <- "
model {
# Likelihood
for(i in 1:n) {
x[i] ~ dnorm(theta, 1/sigma2)
}

# Hierarchical prior
theta ~ dnorm(mu, 1/tau2)
mu ~ dnorm(0, 0.0001)

# Inverse-Gamma prior for tau^2
tau2 ~ dgamma(0.1, 0.1)  
tau <- sqrt(tau2)
}
"

# Prepare data and parameters

data_jags <- list(x = x, n = n, sigma2 = sigma2)
params <- c("theta", "mu", "tau")

# Run JAGS
set.seed(123)
jags_fit <- jags(data = data_jags,
                 parameters.to.save = params,
                 model.file = textConnection(model_string),
                 n.chains = 3,
                 n.iter = 5000,
                 n.burnin = 1000,
                 n.thin = 2)

# Extract posterior chains
posterior_samples <- as.data.frame(jags_fit$BUGSoutput$sims.matrix)
theta_samples <- posterior_samples$theta
mu_samples <- posterior_samples$mu
tau_samples <- posterior_samples$tau

# Compute posterior predictive distribution outside JAGS
n_sims <- length(theta_samples)
x_tilde_samples <- rnorm(n_sims, mean = theta_samples, sd = sigma)

# Posterior summaries
cat("Posterior mean of theta: ", mean(theta_samples), "\n")
cat("Posterior SD of theta:   ", sd(theta_samples), "\n")
cat("Posterior mean of mu:    ", mean(mu_samples), "\n")
cat("Posterior SD of mu:      ", sd(mu_samples), "\n")
cat("Posterior mean of tau:   ", mean(tau_samples), "\n")
cat("Posterior SD of tau:     ", sd(tau_samples), "\n")
cat("Posterior predictive mean of ~X: ", mean(x_tilde_samples), "\n")
cat("Posterior predictive SD of ~X:   ", sd(x_tilde_samples), "\n")

# Plot posterior predictive distribution
hist(x_tilde_samples, breaks=50, col="skyblue",
     main="Posterior Predictive Distribution", xlab=expression(tilde(X)))
abline(v=mean(x_tilde_samples), col="red", lwd=2)
abline(v=mean(x), col="darkgreen", lwd=2, lty=2)
legend("topright", legend=c("Posterior predictive mean","Observed mean"),
       col=c("red","darkgreen"), lwd=2, lty=c(1,2))

