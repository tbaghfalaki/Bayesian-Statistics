# Clear workspace
rm(list = ls())

# Load libraries
library(rstan)
library(sn)  # for the AIS dataset

# Set working directory
setwd("C:/Users/p80744tb/Desktop/hmc")

# Load AIS data
data(ais)

 
# Prepare data for Stan
y <- ais$BMI
x <- as.numeric(ais$sex)-1  # Encode sex: 0=male, 1=female
N <- length(y)

data_list <- list(
  N = N,
  x = x,
  y = y
)

# Set rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Compile Stan model
compiled_model <- stan_model(file = "S1.stan")

# Fit Stan model
fit <- sampling(
  object = compiled_model,
  data = data_list,
  iter = 1000,               # Total iterations per chain
  warmup = 100,              # Burn-in
  chains = 2,                # Number of chains
  cores = 1,                 # Parallel cores
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)

# Print summary of posterior estimates
print(fit, pars = c("beta0", "beta1", "sigma"))

# Diagnostic plots
traceplot(fit, pars = c("beta0", "beta1", "sigma"))
pairs(fit, pars = c("beta0", "beta1", "sigma"))

# Extract posterior samples
posterior_samples <- extract(fit)
head(posterior_samples$beta0)
head(posterior_samples$beta1)
