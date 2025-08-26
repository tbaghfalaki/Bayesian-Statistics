set.seed(123)

# Data
n <- 10
lambda_true <- 2
data <- rexp(n, rate = lambda_true)

# Priors for M1
a <- 0.01
b <- 0.01

# Marginal likelihood for M1
log_marginal_lik_M1 <- a*log(b) - lgamma(a) + lgamma(a + n) - (a + n)*log(b + sum(data))
marginal_lik_M1 <- exp(log_marginal_lik_M1)

# Priors for M2
mu0 <- 0
tau2 <- 10
alpha <- 0.01
beta <- 0.01

# Monte Carlo approximation of marginal likelihood for M2
M <- 1e5
mu_samples <- rnorm(M, mu0, sqrt(tau2))
sigma2_samples <- 1 / rgamma(M, alpha, rate = beta)  # Inv-Gamma sampling

log_lik <- sapply(1:M, function(i) {
  mu <- mu_samples[i]
  sigma2 <- sigma2_samples[i]
  sum(dlnorm(data, meanlog = mu, sdlog = sqrt(sigma2), log = TRUE))
})

marginal_lik_M2 <- mean(exp(log_lik))

# Bayes Factor
BF_12 <- marginal_lik_M1 / marginal_lik_M2

BF_12
log10(BF_12)
