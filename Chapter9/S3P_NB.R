rm(list=ls())
library(rstan)
library(bellreg)
library(dplyr)

setwd("C:\\Users\\p80744tb\\Desktop\\hmc") 


# Load and prepare data
data(cells)

cells <- cells %>%
  mutate(
    smoker = as.numeric(smoker),
    gender = as.numeric(gender),
    weightf = factor(weight, levels = c("normal", "over", "obese")),
    agef = factor(age, levels = c("young", "mid", "old"))
  )

X <- model.matrix(~ smoker + weightf + agef + gender, data = cells)
Y <- cells$cells
N <- nrow(X)
K <- ncol(X)


# Stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# Poisson Model
compiled_model_p <- stan_model(file = "S3P.stan")

stan_data_p <- list(
  N = N, K = K,
  Y = as.integer(Y), X = X  
)

fit_p <- sampling(
  object = compiled_model_p,
  data = stan_data_p,
  iter = 200, warmup = 50, chains = 1, cores = 1,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)

print(fit_p, pars = c("beta"))

# Posterior means
beta_p_mean <- colMeans(extract(fit_p)$beta)
lambda_p_hat <- exp(X %*% beta_p_mean)

# BIC
loglik_p <- sum(dpois(Y, lambda = lambda_p_hat, log = TRUE))
k_p <- length(beta_p_mean)
BIC_p <- log(N) * k_p - 2 * loglik_p
cat("Poisson BIC:", BIC_p, "\n")

# DIC
post_beta_p <- extract(fit_p)$beta
D_p_samples <- sapply(1:nrow(post_beta_p), function(i) {
  -2 * sum(dpois(Y, lambda = exp(X %*% post_beta_p[i,]), log = TRUE))
})
Dbar_p <- mean(D_p_samples)
Dhat_p <- -2 * sum(dpois(Y, lambda = lambda_p_hat, log = TRUE))
DIC_p <- 2 * Dbar_p - Dhat_p
cat("Poisson DIC:", DIC_p, "\n")

# LPML
CPO_p <- sapply(1:N, function(i) {
  mean(1 / dpois(Y[i], lambda = exp(X[i,] %*% t(post_beta_p))))
})
LPML_p <- sum(log(CPO_p))
cat("Poisson LPML:", LPML_p, "\n")


# Negative Binomial Model

compiled_model_nb <- stan_model(file = "S3NB.stan")

stan_data_nb <- list(
  N = N, K = K,
  Y = as.integer(Y), X = X  
)

fit_nb <- sampling(
  object = compiled_model_nb,
  data = stan_data_nb,
  iter = 200, warmup = 50, chains = 1, cores = 1,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)

print(fit_nb, pars = c("beta","phi"))

# Posterior means
posterior_nb <- extract(fit_nb)
beta_nb_mean <- colMeans(posterior_nb$beta)
phi_nb_mean <- mean(posterior_nb$phi)
lambda_nb_hat <- exp(X %*% beta_nb_mean)

# BIC
loglik_nb <- sum(dnbinom(Y, mu = lambda_nb_hat, size = phi_nb_mean, log = TRUE))
k_nb <- length(beta_nb_mean) + 1
BIC_nb <- log(N) * k_nb - 2 * loglik_nb
cat("NB BIC:", BIC_nb, "\n")

# DIC
post_beta_nb <- posterior_nb$beta
post_phi_nb <- posterior_nb$phi
D_nb_samples <- sapply(1:length(post_phi_nb), function(i) {
  -2 * sum(dnbinom(Y, mu = exp(X %*% post_beta_nb[i,]), size = post_phi_nb[i], log = TRUE))
})
Dbar_nb <- mean(D_nb_samples)
Dhat_nb <- -2 * sum(dnbinom(Y, mu = lambda_nb_hat, size = phi_nb_mean, log = TRUE))
DIC_nb <- 2 * Dbar_nb - Dhat_nb
cat("NB DIC:", DIC_nb, "\n")

# LPML
CPO_nb <- sapply(1:N, function(i) {
  mean(1 / dnbinom(Y[i], mu = exp(X[i,] %*% t(post_beta_nb)), size = post_phi_nb))
})
LPML_nb <- sum(log(CPO_nb))
cat("NB LPML:", LPML_nb, "\n")
