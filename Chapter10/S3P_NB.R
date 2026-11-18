
rm(list = ls())

library(rstan)
library(bellreg)
library(MASS)

setwd("/Users/user/Desktop/Lab/Example1")

# Data
data(cells)

cells$smoker <- factor(cells$smoker)
cells$gender <- factor(cells$gender)
cells$weightf <- factor(cells$weight,
                        levels = c("normal", "over", "obese"))
cells$agef <- factor(cells$age,
                     levels = c("young", "mid", "old"))

X <- model.matrix(~ smoker + weightf + agef + gender, data = cells)
Y <- as.integer(cells$cells)
N <- length(Y)
K <- ncol(X)

stan_data <- list(N = N, K = K, Y = Y, X = X)

# Stan settings
rstan_options(auto_write = TRUE)
options(mc.cores = 2)

# Poisson model
fit_p <- stan(
  file = "S3P.stan",
  data = stan_data,
  chains = 2, iter = 2000, warmup = 1000,
  seed = 1234
)

print(fit_p, pars = "beta")

# Negative Binomial model
fit_nb <- stan(
  file = "S3NB.stan",
  data = stan_data,
  chains = 2, iter = 2000, warmup = 1000,
  seed = 1234 
)

print(fit_nb, pars = c("beta", "phi"))

# Posterior samples
post_p <- extract(fit_p)
post_nb <- extract(fit_nb)

beta_p <- post_p$beta
beta_nb <- post_nb$beta
phi_nb <- post_nb$phi

# ------------------------------------------------
# BIC 
# ------------------------------------------------

# Poisson
loglik_p <- max(rowSums(post_p$log_lik))
BIC_p <- -2 * loglik_p + K * log(N)

# Negative Binomial
loglik_nb <- max(rowSums(post_nb$log_lik))
BIC_nb <- -2 * loglik_nb + (K + 1) * log(N)

cat("Poisson BIC:", BIC_p, "\n")
cat("Negative Binomial BIC:", BIC_nb, "\n")


# ------------------------------------------------
# DIC
# ------------------------------------------------

# Poisson
Dbar_p <- mean(-2 * rowSums(post_p$log_lik))

mu_p <- exp(X %*% colMeans(beta_p))
Dhat_p <- -2 * sum(dpois(Y, mu_p, log = TRUE))

DIC_p <- 2 * Dbar_p - Dhat_p

# Negative Binomial
Dbar_nb <- mean(-2 * rowSums(post_nb$log_lik))

mu_nb <- exp(X %*% colMeans(beta_nb))
phi_hat <- mean(phi_nb)

Dhat_nb <- -2 * sum(dnbinom(Y, mu = mu_nb,
                            size = phi_hat, log = TRUE))

DIC_nb <- 2 * Dbar_nb - Dhat_nb

# ------------------------------------------------
# LPML
# ------------------------------------------------

# Poisson
LPML_p <- -sum(log(colMeans(exp(-post_p$log_lik))))

# Negative Binomial
LPML_nb <- -sum(log(colMeans(exp(-post_nb$log_lik))))

# ------------------------------------------------
# Model comparison
# ------------------------------------------------

results <- data.frame(
  Model = c("Poisson", "Negative Binomial"),
  BIC = c(BIC_p, BIC_nb),
  DIC = c(DIC_p, DIC_nb),
  LPML = c(LPML_p, LPML_nb)
)

print(results)



library(bayesplot)

# Trace plots
mcmc_trace(as.array(fit_p), regex_pars = "^beta\\[")

# Posterior density plots
mcmc_dens_overlay(as.array(fit_p), regex_pars = "^beta\\[")

# Autocorrelation plots
mcmc_acf(as.array(fit_p), regex_pars = "^beta\\[", lags = 30)




# Trace plots
mcmc_trace(as.array(fit_nb), regex_pars = "^(beta\\[|phi$)")

# Posterior density plots
mcmc_dens_overlay(as.array(fit_nb), regex_pars = "^(beta\\[|phi$)")

# Autocorrelation plots
mcmc_acf(as.array(fit_nb), regex_pars = "^(beta\\[|phi$)",
         lags = 30)
