rm(list = ls())

library(ggplot2)
library(survminer)
library(rstan)
library(JM)

# Load AIDS dataset
data(aids.id)

setwd("C:\\Users\\p80744tb\\Desktop\\hmc\\survival") 

# Convert categorical covariates to numeric (0/1)
aids.id$drug   <- ifelse(aids.id$drug == "ddI", 1, 0)
aids.id$gender <- ifelse(aids.id$gender == "male", 1, 0)
aids.id$prevOI <- ifelse(aids.id$prevOI == "AIDS", 1, 0) 
aids.id$AZT    <- ifelse(aids.id$AZT == "intolerance", 1, 0)

# Design matrix (including intercept)
X <- model.matrix(~ drug + gender + prevOI + AZT, data = aids.id)
N <- nrow(aids.id)
K <- ncol(X)

# Piecewise-constant baseline hazard: define intervals
breaks <- seq(0, 18, by = 3)  # e.g., 0,3,6,9,12,15,18
J <- length(breaks) - 1       # number of intervals

# Compute interval index for each observation
get_interval <- function(t) {
  i <- findInterval(t, breaks, rightmost.closed = TRUE)
  if (i > J) i <- J           # cap at last interval
  return(i)
}
interval <- sapply(aids.id$Time, get_interval)

# Prepare data for Stan (include breaks)
stan_data <- list(
  N = N,
  K = K,
  J = J,
  X = X,
  t = aids.id$Time,
  delta = aids.id$death,
  interval = interval,
  breaks = breaks              # <--- must include this for Stan
)

# Set rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Compile Stan model
compiled_model <- stan_model(file = "CoxPHpw.stan")

# Sampling with Stan
fit_PH <- sampling(
  object = compiled_model,
  data = stan_data,
  chains = 2,
  iter = 1000,
  warmup = 500,
  seed = 123
)

# Print posterior summaries
print(fit_PH, pars = c("beta","lambda"))



# Extract log-likelihood per draw
log_lik <- extract(fit_PH, "log_lik")$log_lik  # draws x N

# Posterior mean log-likelihood per observation
mean_loglik <- apply(log_lik, 2, mean)

# DIC
D_bar <- -2 * sum(mean_loglik)                     # deviance at posterior mean
D_hat <- -2 * sum(log(colMeans(exp(log_lik))))     # expected deviance
pD <- D_bar - D_hat                                 # effective # parameters
DIC <- D_bar + pD
DIC

# BIC
n <- nrow(X)
p <- length(fit_PH@model_pars)            # number of parameters
logL_hat <- sum(mean_loglik)
BIC <- -2 * logL_hat + p * log(n)
BIC

# LPML (via harmonic mean of conditional predictive ordinates)
CPO <- 1 / apply(exp(-log_lik), 2, mean)
LPML <- sum(log(CPO))
LPML

