
rm(list = ls())

library(rstan)
library(bayesplot)

setwd("/Users/user/Desktop/Lab/Example2")

# ------------------------------------------------
# Data
# ------------------------------------------------

Data <- read.table("BHPS.txt", header = TRUE)

Data$age <- as.numeric(scale(Data$age))

y <- Data$income.1

X <- model.matrix(
  ~ sex + age + mar1 + mar2 + mar3 +
    edu1 + edu2 + edu3 + leisure.1,
  data = Data
)

N <- length(y)
K <- ncol(X)

stopifnot(nrow(X) == N, all(is.finite(y)), all(y > 0))

stan_data <- list(N = N, K = K, y = y, X = X)

cat("N =", N, "; K =", K, "\n")
print(colnames(X))

# ------------------------------------------------
# Stan settings
# ------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = 2)

# ------------------------------------------------
# Fit the three models
# ------------------------------------------------

# Exponential
fit_exp <- stan(
  file = "EXP.stan",
  data = stan_data,
  chains = 2,
  iter = 2000,
  warmup = 1000,
  seed = 1234,
  cores = 2,
  control = list(adapt_delta = 0.95)
)

# Weibull
fit_weib <- stan(
  file = "WEIBULL.stan",
  data = stan_data,
  chains = 2,
  iter = 2000,
  warmup = 1000,
  seed = 2345,
  cores = 2,
  control = list(adapt_delta = 0.95)
)

# Lognormal
fit_logn <- stan(
  file = "LOGNORMAL.stan",
  data = stan_data,
  chains = 2,
  iter = 2000,
  warmup = 1000,
  seed = 3456,
  cores = 2,
  control = list(adapt_delta = 0.95)
)

# Posterior summaries
print(fit_exp, pars = "beta")
print(fit_weib, pars = c("beta", "shape"))
print(fit_logn, pars = c("beta", "sigma"))

# Save fitted models
saveRDS(fit_exp, "fit_exp.rds")
saveRDS(fit_weib, "fit_weib.rds")
saveRDS(fit_logn, "fit_logn.rds")

# ------------------------------------------------
# Posterior samples
# ------------------------------------------------

post_exp <- extract(fit_exp)
post_weib <- extract(fit_weib)
post_logn <- extract(fit_logn)

# Calculate log-likelihood from posterior draws
# if it was not generated in Stan.

# Exponential
if (is.null(post_exp$log_lik)) {
  post_exp$log_lik <- t(sapply(
    1:nrow(post_exp$beta), function(s) {
      eta <- drop(X %*% post_exp$beta[s, ])
      dexp(y, rate = exp(-eta), log = TRUE)
    }
  ))
}

# Weibull
if (is.null(post_weib$log_lik)) {
  post_weib$log_lik <- t(sapply(
    1:nrow(post_weib$beta), function(s) {
      eta <- drop(X %*% post_weib$beta[s, ])
      dweibull(y, shape = post_weib$shape[s],
               scale = exp(eta), log = TRUE)
    }
  ))
}

# Lognormal
if (is.null(post_logn$log_lik)) {
  post_logn$log_lik <- t(sapply(
    1:nrow(post_logn$beta), function(s) {
      eta <- drop(X %*% post_logn$beta[s, ])
      dlnorm(y, meanlog = eta,
             sdlog = post_logn$sigma[s], log = TRUE)
    }
  ))
}

# Check dimensions
print(dim(post_exp$log_lik))
print(dim(post_weib$log_lik))
print(dim(post_logn$log_lik))

# ------------------------------------------------
# BIC (approximation using Stan draws)
# ------------------------------------------------

BIC_exp <- -2 * max(rowSums(post_exp$log_lik)) +
  K * log(N)

BIC_weib <- -2 * max(rowSums(post_weib$log_lik)) +
  (K + 1) * log(N)

BIC_logn <- -2 * max(rowSums(post_logn$log_lik)) +
  (K + 1) * log(N)

# ------------------------------------------------
# DIC
# ------------------------------------------------

# Posterior mean deviance
Dbar_exp <- mean(-2 * rowSums(post_exp$log_lik))
Dbar_weib <- mean(-2 * rowSums(post_weib$log_lik))
Dbar_logn <- mean(-2 * rowSums(post_logn$log_lik))

# Linear predictors at posterior means
eta_exp <- drop(X %*% colMeans(post_exp$beta))
eta_weib <- drop(X %*% colMeans(post_weib$beta))
eta_logn <- drop(X %*% colMeans(post_logn$beta))

# Deviance at posterior mean
Dhat_exp <- -2 * sum(
  dexp(y, rate = exp(-eta_exp), log = TRUE)
)

Dhat_weib <- -2 * sum(
  dweibull(y,
           shape = mean(post_weib$shape),
           scale = exp(eta_weib),
           log = TRUE)
)

Dhat_logn <- -2 * sum(
  dlnorm(y,
         meanlog = eta_logn,
         sdlog = mean(post_logn$sigma),
         log = TRUE)
)

# DIC
DIC_exp <- 2 * Dbar_exp - Dhat_exp
DIC_weib <- 2 * Dbar_weib - Dhat_weib
DIC_logn <- 2 * Dbar_logn - Dhat_logn

# ------------------------------------------------
# LPML
# ------------------------------------------------

# Numerically stable log-mean-exp
log_mean_exp <- function(x) {
  a <- max(x)
  a + log(mean(exp(x - a)))
}

LPML_exp <- -sum(
  apply(-post_exp$log_lik, 2, log_mean_exp)
)

LPML_weib <- -sum(
  apply(-post_weib$log_lik, 2, log_mean_exp)
)

LPML_logn <- -sum(
  apply(-post_logn$log_lik, 2, log_mean_exp)
)

# ------------------------------------------------
# Model comparison
# ------------------------------------------------

results <- data.frame(
  Model = c("Exponential", "Weibull", "Lognormal"),
  BIC = c(BIC_exp, BIC_weib, BIC_logn),
  DIC = c(DIC_exp, DIC_weib, DIC_logn),
  LPML = c(LPML_exp, LPML_weib, LPML_logn)
)

print(results, row.names = FALSE, digits = 6)

# Smaller BIC and DIC are preferred.
# Larger LPML is preferred.

# ------------------------------------------------
# MCMC convergence diagnostics
# ------------------------------------------------

print(summary(fit_exp)$summary[, c("n_eff", "Rhat")])
print(summary(fit_weib)$summary[, c("n_eff", "Rhat")])
print(summary(fit_logn)$summary[, c("n_eff", "Rhat")])

# Divergent transitions
for (fit in list(fit_exp, fit_weib, fit_logn)) {
  sp <- get_sampler_params(fit, inc_warmup = FALSE)
  
  n_div <- sum(sapply(
    sp, function(s) sum(s[, "divergent__"])
  ))
  
  cat("Divergences:", n_div, "\n")
}

# ------------------------------------------------
# Diagnostic plots
# ------------------------------------------------

color_scheme_set("blue")

# Exponential
draws_exp <- as.array(fit_exp)

print(mcmc_trace(draws_exp, regex_pars = "^beta\\["))
print(mcmc_dens_overlay(draws_exp, regex_pars = "^beta\\["))
print(mcmc_acf(draws_exp, regex_pars = "^beta\\[",
               lags = 30))

# Weibull
draws_weib <- as.array(fit_weib)

print(mcmc_trace(draws_weib,
                 regex_pars = "^(beta\\[|shape$)"))

print(mcmc_dens_overlay(draws_weib,
                        regex_pars = "^(beta\\[|shape$)"))

print(mcmc_acf(draws_weib,
               regex_pars = "^(beta\\[|shape$)",
               lags = 30))

# Lognormal
draws_logn <- as.array(fit_logn)

print(mcmc_trace(draws_logn,
                 regex_pars = "^(beta\\[|sigma$)"))

print(mcmc_dens_overlay(draws_logn,
                        regex_pars = "^(beta\\[|sigma$)"))

print(mcmc_acf(draws_logn,
               regex_pars = "^(beta\\[|sigma$)",
               lags = 30))

