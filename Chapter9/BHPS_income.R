rm(list = ls())

library(rstan)
library(dplyr)

setwd("C:\\Users\\p80744tb\\Desktop\\hmc\\glm")

# ---- Load Data ----
Data <- read.table("BHPS.txt", header = TRUE)

# ---- Prepare data ----
y <- Data$income.1
Data$age=(Data$age-mean(Data$age))/sd(Data$age)
X <- model.matrix(
  ~ age + sex + edu1 + edu2 + edu3 + mar1 + mar2 + mar3 + leisure.1,
  data = Data
)

N <- nrow(X)
K <- ncol(X)

# Stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ---- Compile model ----
compiled_exp <- stan_model(file = "EXP.stan")

stan_data_exp <- list(
  N = N,
  K = K,
  y = y,
  X = X
)
#lm(income.1  ~ age + sex + edu1 + edu2 + edu3 + mar1 + mar2 + mar3 + leisure.1, data = Data)
# ---- Fit model ----
fit_exp <- sampling(
  object = compiled_exp,
  data = stan_data_exp,
  iter = 200, warmup = 100,
  chains = 1, cores = 1,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

print(fit_exp, pars = "beta")

# ---- Posterior means ----
posterior_exp <- extract(fit_exp)
beta_exp_mean <- colMeans(posterior_exp$beta)

eta_hat <- X %*% beta_exp_mean
rate_hat <- exp(eta_hat)

# ----------------------------
#   BIC
# ----------------------------
loglik_exp <- sum(dexp(y, rate = rate_hat, log = TRUE))

k_exp <- length(beta_exp_mean)
BIC_exp <- log(N) * k_exp - 2 * loglik_exp
cat("Exponential BIC:", BIC_exp, "\n")

# ----------------------------
#   DIC
# ----------------------------
post_beta <- posterior_exp$beta

D_exp_samples <- sapply(1:nrow(post_beta), function(i) {
  rate_i <- exp(X %*% post_beta[i,])
  -2 * sum(dexp(y, rate = rate_i, log = TRUE))
})

Dbar_exp <- mean(D_exp_samples)
Dhat_exp <- -2 * sum(dexp(y, rate = rate_hat, log = TRUE))

DIC_exp <- 2 * Dbar_exp - Dhat_exp
cat("Exponential DIC:", DIC_exp, "\n")

# ----------------------------
#   LPML
# ----------------------------
CPO_exp <- sapply(1:N, function(i) {
  mean(1 / dexp(y[i], rate = exp(X[i,] %*% t(post_beta))))
})

LPML_exp <- sum(log(CPO_exp))
cat("Exponential LPML:", LPML_exp, "\n")

###############
#################
# ---- Compile model ----
compiled_gamma <- stan_model(file = "GAMMA.stan")

stan_data_gamma <- list(
  N = N,
  K = K,
  y = y,
  X = X
)

# ---- Fit model ----
fit_gamma <- sampling(
  object = compiled_gamma,
  data = stan_data_gamma,
  iter = 200, warmup = 100,
  chains = 1, cores = 1,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

print(fit_gamma, pars = c("beta", "alpha"))

# ---- Posterior means ----
posterior_gamma <- extract(fit_gamma)
beta_gamma_mean <- colMeans(posterior_gamma$beta)
tau_gamma_mean <- mean(posterior_gamma$alpha)

mu_hat_gamma <- exp(X %*% beta_gamma_mean)

# ----------------------------
#   BIC
# ----------------------------
# Gamma parameterization: scale = mu / tau
loglik_gamma <- sum(dgamma(y, shape = mu_hat_gamma * tau_gamma_mean, scale = 1 / tau_gamma_mean, log = TRUE))
k_gamma <- length(beta_gamma_mean) + 1  # +1 for tau
BIC_gamma <- log(N) * k_gamma - 2 * loglik_gamma
cat("Gamma BIC:", BIC_gamma, "\n")

# ----------------------------
#   DIC
# ----------------------------
post_beta <- posterior_gamma$beta
post_alpha <- posterior_gamma$alpha

D_gamma_samples <- sapply(1:length(post_alpha), function(i) {
  mu_i <- exp(X %*% post_beta[i,])
  -2 * sum(dgamma(y, shape = mu_i * post_alpha[i], scale = 1 / post_alpha[i], log = TRUE))
})

Dbar_gamma <- mean(D_gamma_samples)
Dhat_gamma <- -2 * sum(dgamma(y, shape = mu_hat_gamma * tau_gamma_mean, scale = 1 / tau_gamma_mean, log = TRUE))
DIC_gamma <- 2 * Dbar_gamma - Dhat_gamma
cat("Gamma DIC:", DIC_gamma, "\n")

# ----------------------------
#   LPML
# ----------------------------
 
CPO_gamma <- sapply(1:N, function(i) {
  eta_i <- as.numeric(X[i, ] %*% t(post_beta))  # 1 x S vector
  shape_i <- eta_i * post_alpha                    # element-wise multiplication
  mean(1 / dgamma(y[i], shape = shape_i, scale = 1 / post_alpha))
})



LPML_gamma <- sum(log(CPO_gamma))
cat("Gamma LPML:", LPML_gamma, "\n")






#################
# ---- Compile model ----
compiled_weib <- stan_model(file = "WEIBULL.stan")

stan_data_weib <- list(
  N = N,
  K = K,
  y = y,
  X = X
)

# ---- Fit model ----
fit_weib <- sampling(
  object = compiled_weib,
  data = stan_data_weib,
  iter = 200, warmup = 100,
  chains = 1, cores = 1,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

print(fit_weib, pars = c("beta", "shape"))

# ---- Posterior means ----
posterior_weib <- extract(fit_weib)
beta_weib_mean <- colMeans(posterior_weib$beta)
shape_weib_mean <- mean(posterior_weib$shape)

scale_hat <- exp(X %*% beta_weib_mean)

# ----------------------------
#   BIC
# ----------------------------
loglik_weib <- sum(dweibull(y, shape = shape_weib_mean, scale = scale_hat, log = TRUE))
k_weib <- length(beta_weib_mean) + 1  # +1 for shape
BIC_weib <- log(N) * k_weib - 2 * loglik_weib
cat("Weibull BIC:", BIC_weib, "\n")

# ----------------------------
#   DIC
# ----------------------------
post_beta <- posterior_weib$beta
post_shape <- posterior_weib$shape

D_weib_samples <- sapply(1:length(post_shape), function(i) {
  scale_i <- exp(X %*% post_beta[i,])
  -2 * sum(dweibull(y, shape = post_shape[i], scale = scale_i, log = TRUE))
})

Dbar_weib <- mean(D_weib_samples)
Dhat_weib <- -2 * sum(dweibull(y, shape = shape_weib_mean, scale = scale_hat, log = TRUE))
DIC_weib <- 2 * Dbar_weib - Dhat_weib
cat("Weibull DIC:", DIC_weib, "\n")

# ----------------------------
#   LPML
# ----------------------------
CPO_weib <- sapply(1:N, function(i) {
  mean(1 / dweibull(y[i], shape = post_shape, scale = exp(X[i,] %*% t(post_beta))))
})

LPML_weib <- sum(log(CPO_weib))
cat("Weibull LPML:", LPML_weib, "\n")


#########################
# Stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ---- Compile model ----
compiled_logn <- stan_model(file = "LOGNORMAL.stan")

stan_data_logn <- list(
  N = N,
  K = K,
  y = y,
  X = X
)

# ---- Fit model ----
fit_logn <- sampling(
  object = compiled_logn,
  data = stan_data_logn,
  iter = 200, warmup = 100,
  chains = 1, cores = 1,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

print(fit_logn, pars = c("beta", "sigma"))

# ---- Posterior means ----
posterior_logn <- extract(fit_logn)
beta_logn_mean <- colMeans(posterior_logn$beta)
sigma_logn_mean <- mean(posterior_logn$sigma)

mu_hat <- X %*% beta_logn_mean

# ----------------------------
#   BIC
# ----------------------------
loglik_logn <- sum(dlnorm(y, meanlog = mu_hat, sdlog = sigma_logn_mean, log = TRUE))
k_logn <- length(beta_logn_mean) + 1  # +1 for sigma
BIC_logn <- log(N) * k_logn - 2 * loglik_logn
cat("Lognormal BIC:", BIC_logn, "\n")

# ----------------------------
#   DIC
# ----------------------------
post_beta <- posterior_logn$beta
post_sigma <- posterior_logn$sigma

D_logn_samples <- sapply(1:length(post_sigma), function(i) {
  -2 * sum(dlnorm(y, meanlog = X %*% post_beta[i,], sdlog = post_sigma[i], log = TRUE))
})

Dbar_logn <- mean(D_logn_samples)
Dhat_logn <- -2 * sum(dlnorm(y, meanlog = mu_hat, sdlog = sigma_logn_mean, log = TRUE))
DIC_logn <- 2 * Dbar_logn - Dhat_logn
cat("Lognormal DIC:", DIC_logn, "\n")

# ----------------------------
#   LPML
# ----------------------------
CPO_logn <- sapply(1:N, function(i) {
  mean(1 / dlnorm(y[i], meanlog = X[i,] %*% t(post_beta), sdlog = post_sigma))
})

LPML_logn <- sum(log(CPO_logn))
cat("Lognormal LPML:", LPML_logn, "\n")









print(fit_exp, pars = "beta")
cat("Exponential BIC:", BIC_exp, "\n")
cat("Exponential DIC:", DIC_exp, "\n")
cat("Exponential LPML:", LPML_exp, "\n")


print(fit_logn, pars = c("beta", "sigma"))
cat("Lognormal BIC:", BIC_logn, "\n")
cat("Lognormal DIC:", DIC_logn, "\n")
cat("Lognormal LPML:", LPML_logn, "\n")


print(fit_weib, pars = c("beta", "shape"))
cat("Weibull BIC:", BIC_weib, "\n")
cat("Weibull DIC:", DIC_weib, "\n")
cat("Weibull LPML:", LPML_weib, "\n")





print(fit_gamma, pars = c("beta", "tau"))
cat("gamma BIC:", BIC_gamma, "\n")
cat("Gamma DIC:", DIC_gamma, "\n")
cat("Gamma LPML:", LPML_gamma, "\n")