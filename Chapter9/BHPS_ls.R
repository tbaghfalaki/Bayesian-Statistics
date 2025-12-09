rm(list=ls())
# Load required libraries
library(rstan)
library(dplyr)

# Read the data
Data <- read.table("BHPS.txt", header = TRUE)

# Select variables for the ordinal model
y <- Data$lifesatisfaction.1
X <- Data %>%
  select(age, sex, edu1, edu2, edu3, mar1, mar2, mar3, leisure.1) %>%
  as.matrix()

# Number of covariates
P <- ncol(X)
# Number of observations
N <- nrow(X)
# Number of ordinal categories
K <- length(unique(y))

# Prepare data for Stan
stan_data <- list(
  N = N,
  K = K,
  P = P,
  X = X,
  y = y
)

# Stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ---- Compile model ----
compiled_ordinal <- stan_model(file = "ordinal.stan")

# ---- Fit model ----
fit_ordinal <- sampling(
  object = compiled_ordinal,
  data = stan_data,
  iter = 200, warmup = 100,
  chains = 1, cores = 1,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  seed = 123
)

# ---- Print posterior summaries ----
print(fit_ordinal, pars = c("beta", "c"), probs = c(0.025, 0.5, 0.975))

