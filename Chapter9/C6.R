rm(list = ls())

library(ggplot2)
library(survminer)
library(rstan)
library(JM)

# Load data
data(aids.id)
names(aids.id)

# Set working directory (adjust path)
setwd("C:\\Users\\p80744tb\\Desktop\\hmc\\survival") 

# Convert categorical variables to numeric indicators
aids.id$drug   <- ifelse(aids.id$drug == "ddI", 1, 0)
aids.id$gender <- ifelse(aids.id$gender == "male", 1, 0)
aids.id$prevOI <- ifelse(aids.id$prevOI == "AIDS", 1, 0) 
aids.id$AZT    <- ifelse(aids.id$AZT == "intolerance", 1, 0)

# Design matrix (no intercept; lambda0 absorbs intercept)
X <- as.matrix(aids.id[, c("drug","gender","prevOI","AZT")])
N <- nrow(X)
K <- ncol(X)

# Stan data
stan_data <- list(
  N = N,
  K = K,
  X = X,
  T = aids.id$Time,
  delta = aids.id$death
)

# rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Compile Stan model
compiled_model <- stan_model(file = "CoxPH.stan")

# Sampling
fit_PH <- sampling(
  object = compiled_model,
  data = stan_data,
  chains = 1,
  iter = 400,
  warmup = 200,
  thin = 1,
  seed = 1234,
  cores = 4
)

# Print posterior summaries
print(fit_PH, pars = c("beta","lambda0"))

# Extract log-likelihood
log_lik <- extract(fit_PH, "log_lik")$log_lik  # draws x N

# Posterior mean log-likelihood per observation
mean_loglik <- apply(log_lik, 2, mean)

# DIC
D_bar <- -2 * sum(mean_loglik)
D_hat <- -2 * sum(log(colMeans(exp(log_lik))))
pD <- D_bar - D_hat
DIC4 <- D_bar + pD
DIC4

# BIC
logL_hat <- sum(mean_loglik)
p <- K + 1    # K coefficients + lambda0
BIC4 <- -2 * logL_hat + p * log(N)
BIC4

# LPML via harmonic mean approximation
CPO <- 1 / apply(exp(-log_lik), 2, mean)
LPML4 <- sum(log(CPO))
LPML4
#######################