rm(list=ls())
library(rstan)
library(bellreg)
library(dplyr)

setwd("C:\\Users\\p80744tb\\Desktop\\hmc") 

#---------------------------
# Load and prepare data
#---------------------------
data(cells)

cells <- cells %>%
  mutate(
    smoker = as.numeric(smoker),
    gender = as.numeric(gender),
    weightf = factor(weight, levels = c("normal", "over", "obese")),
    agef = factor(age, levels = c("young", "mid", "old"))
  )

X <- model.matrix(~ smoker + weightf + agef + gender, data = cells)
Z <- X
Y <- cells$cells
N <- nrow(X)
K <- ncol(X)
Q <- ncol(Z)

#---------------------------
# Helper function for LPML
#---------------------------
compute_lpml <- function(log_lik_mat) {
  max_loglik <- apply(log_lik_mat, 2, max)
  lik_mat <- exp(log_lik_mat - matrix(max_loglik, nrow = nrow(log_lik_mat),
                                      ncol = ncol(log_lik_mat), byrow = TRUE))
  
  cpo <- numeric(N)
  for(i in 1:N) {
    cpo[i] <- 1 / mean(1 / lik_mat[, i])
  }
  lpml <- sum(log(cpo)) + sum(max_loglik)
  return(lpml)
}

#---------------------------
# Stan options
#---------------------------
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#---------------------------
# ZIP Model
#---------------------------
stan_data_zip <- list(
  N = N, K = K, Q = Q,
  Y = as.integer(Y), X = X, Z = Z,
  sigma_beta = 5, sigma_gamma = 5
)

compiled_model_zip <- stan_model(file = "S3ZIP.stan")

fit_zip <- sampling(
  object = compiled_model_zip,
  data = stan_data_zip,
  iter = 200, warmup = 50, chains = 1, cores = 1,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)

print(fit_zip, pars = c("beta", "gamma"))

log_lik_zip <- extract(fit_zip, pars = "log_lik")$log_lik

# DIC
D_zip <- -2 * log_lik_zip
D_bar_zip <- mean(rowSums(D_zip))
D_hat_zip <- -2 * sum(log_lik_zip[1,])
pD_zip <- D_bar_zip - D_hat_zip
DIC_zip <- D_bar_zip + pD_zip

# BIC
beta_mean_zip <- colMeans(extract(fit_zip, "beta")$beta)
gamma_mean_zip <- colMeans(extract(fit_zip, "gamma")$gamma)
loglik_at_mean_zip <- sum(dpois(Y, lambda = exp(X %*% beta_mean_zip), log = TRUE))
k_zip <- length(beta_mean_zip) + length(gamma_mean_zip)
BIC_zip <- -2 * loglik_at_mean_zip + k_zip * log(N)

# LPML
LPML_zip <- compute_lpml(log_lik_zip)

#---------------------------
# ZINB Model
#---------------------------
stan_data_zinb <- list(
  N = N, K = K, Q = Q,
  Y = as.integer(Y), X = X, Z = Z,
  sigma_beta = 5, sigma_gamma = 5,
  a_phi = 2, b_phi = 0.1
)

compiled_model_zinb <- stan_model(file = "S3ZINB.stan")

fit_zinb <- sampling(
  object = compiled_model_zinb,
  data = stan_data_zinb,
  iter = 200, warmup = 50, chains = 1, cores = 1,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)

print(fit_zinb, pars = c("beta", "gamma", "phi"))

log_lik_zinb <- extract(fit_zinb, pars = "log_lik")$log_lik

# DIC
D_zinb <- -2 * log_lik_zinb
D_bar_zinb <- mean(rowSums(D_zinb))
D_hat_zinb <- -2 * sum(log_lik_zinb[1,])
pD_zinb <- D_bar_zinb - D_hat_zinb
DIC_zinb <- D_bar_zinb + pD_zinb

# BIC
beta_mean_zinb <- colMeans(extract(fit_zinb, "beta")$beta)
gamma_mean_zinb <- colMeans(extract(fit_zinb, "gamma")$gamma)
loglik_at_mean_zinb <- sum(dpois(Y, lambda = exp(X %*% beta_mean_zinb), log = TRUE))
k_zinb <- length(beta_mean_zinb) + length(gamma_mean_zinb)
BIC_zinb <- -2 * loglik_at_mean_zinb + k_zinb * log(N)

# LPML
LPML_zinb <- compute_lpml(log_lik_zinb)

#---------------------------
# Print results
#---------------------------
cat("ZIP Model:\nDIC:", DIC_zip, "\nBIC:", BIC_zip, "\nLPML:", LPML_zip, "\n\n")
cat("ZINB Model:\nDIC:", DIC_zinb, "\nBIC:", BIC_zinb, "\nLPML:", LPML_zinb, "\n")
