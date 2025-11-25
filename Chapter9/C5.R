rm(list=ls())
library(ggplot2)
library(survminer)
library(rstan)
library(JM)

data(aids.id)
names(aids.id)

setwd("C:\\Users\\p80744tb\\Desktop\\hmc\\survival") 

Bfit1 <- survfit(Surv(Time, death) ~ drug, data = aids.id) 
ggsurvplot(Bfit1, conf.int = TRUE,
           risk.table = FALSE, risk.table.col = "strata",
           ggtheme = theme_bw() , legend.title = element_blank(), 
           legend.labs = c("ddC", "ddI"),
           palette = c("#1F77B4", "#FF7F0E") ) 





 
# Convert to numeric indicators
aids.id$drug   <- ifelse(aids.id$drug == "ddI", 1, 0)
aids.id$gender <- ifelse(aids.id$gender == "male", 1, 0)
aids.id$prevOI <- ifelse(aids.id$prevOI == "AIDS", 1, 0) 
aids.id$AZT    <- ifelse(aids.id$AZT == "intolerance", 1, 0)

# Log-time
logT <- log(aids.id$Time)

# Design matrix
X <- model.matrix(~ drug + gender + prevOI + AZT, data = aids.id)

stan_data <- list(
  N = nrow(aids.id),
  logT = logT,
  delta = aids.id$death,
  K = ncol(X),
  X = X
)

# Set rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# Compile Stan model
compiled_model <- stan_model(file = "LN.stan")


# Sampling with Stan
fit_LN <- sampling(
  object = compiled_model,
  data = stan_data,
  chains = 1,
  iter = 400,
  warmup = 200,
  thin = 1,
  seed = 1234,
  cores = 4
)


print(fit_LN, pars = c("beta","sigma"))

# Extract log-likelihood per draw
log_lik <- extract(fit_LN, "log_lik")$log_lik  # draws x N

# Posterior mean log-likelihood per observation
mean_loglik <- apply(log_lik, 2, mean)

# DIC
D_bar <- -2 * sum(mean_loglik)                     # deviance at posterior mean
D_hat <- -2 * sum(log(colMeans(exp(log_lik))))     # expected deviance
pD <- D_bar - D_hat                                 # effective # parameters
DIC1 <- D_bar + pD
DIC1

# BIC
n <- ncol(X)
p <- length(fit_LN@model_pars)            # number of parameters
logL_hat <- sum(mean_loglik)
BIC1 <- -2 * logL_hat + p * log(n)
BIC1

# LPML via harmonic mean approximation
CPO <- 1 / apply(exp(-log_lik), 2, mean)
LPML1 <- sum(log(CPO))
LPML1



###########################################################################
stan_data <- list(
  N = nrow(aids.id),
  T = aids.id$Time,
  delta = aids.id$death,
  K = ncol(X),
  X = X
)
# Set rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# Compile Stan model
compiled_model <- stan_model(file = "LL.stan")


# Sampling with Stan
fit_LL <- sampling(
  object = compiled_model,
  data = stan_data,
  chains = 1,
  iter = 400,
  warmup = 200,
  thin = 1,
  seed = 1234,
  cores = 4
)


print(fit_LL, pars = c("beta","sigma"))


# Extract log-likelihood per draw
log_lik <- extract(fit_LL, "log_lik")$log_lik  # dim: draws x N

# Posterior mean log-likelihood per obs
mean_loglik <- apply(log_lik, 2, mean)

# --- DIC ---
D_bar <- -2 * sum(mean_loglik)             # mean deviance
D_hat <- -2 * sum(log(colMeans(exp(log_lik))))  # deviance at posterior mean
pD <- D_bar - D_hat                         # effective # parameters
DIC2 <- D_bar + pD
DIC2

# --- BIC ---
n <- ncol(X)
p <- length(fit_LL@model_pars)  # number of parameters
logL_hat <- sum(mean_loglik)
BIC2 <- -2 * logL_hat + p * log(n)
BIC2

# --- LPML (using harmonic mean approx for CPO) ---
CPO <- 1 / apply(exp(-log_lik), 2, mean)
LPML2 <- sum(log(CPO))
LPML2

###########################################################################

stan_data <- list(
  N = nrow(aids.id),
  T = aids.id$Time,
  delta = aids.id$death,
  K = ncol(X),
  X = X
)
# Set rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# Compile Stan model
compiled_model <- stan_model(file = "Weibull.stan")


# Sampling with Stan
fit_weibull <- sampling(
  object = compiled_model,
  data = stan_data,
  chains = 1,
  iter = 400,
  warmup = 200,
  thin = 1,
  seed = 1234,
  cores = 4
)


print(fit_weibull, pars = c("beta","alpha"))


# Extract log-likelihood per draw
log_lik <- extract(fit_weibull, "log_lik")$log_lik  # dim: draws x N

# Posterior mean log-likelihood per obs
mean_loglik <- apply(log_lik, 2, mean)

# DIC
D_hat <- -2 * sum(log(mean(exp(log_lik)) ))        # expected deviance
D_bar <- -2 * sum(mean_loglik)                     # deviance at posterior mean
pD <- D_bar - D_hat                                 # effective number of parameters
DIC3 <- D_bar + pD
DIC3

# BIC approximation
p <- length(fit_weibull@model_pars)   # number of parameters
logL_hat <- sum(mean_loglik)
BIC3 <- -2*logL_hat + p*log(nrow(aids.id))
BIC3

# LPML via CPO (using harmonic mean approximation)
CPO <- 1 / apply(exp(-log_lik), 2, mean)
LPML3 <- sum(log(CPO))
LPML3

###########################################################################
DIC1
DIC2
DIC3

BIC1
BIC2
BIC3

LPML1
LPML2
LPML3
