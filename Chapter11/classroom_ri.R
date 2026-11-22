library(rstan)
library(dplyr)
library(nlme)
setwd('C:\\Users\\p80744tb\\Desktop\\code')

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# load and prepare data
data(MathAchieve)
df <- as.data.frame(MathAchieve)

# create numeric indicators
df$sex      <- ifelse(df$Sex      == "Female", 1, 0)
df$minority <- ifelse(df$Minority == "Yes",    1, 0)
df$school   <- as.integer(factor(df$School))

# fixed-effects design matrix: intercept + SES + minority + sex
X <- model.matrix(~ SES + minority + sex, data = df)

stan_data <- list(
  N      = nrow(df),
  J      = length(unique(df$school)),
  K      = ncol(X),
  X      = X,
  school = df$school,
  y      = df$MathAch
)

fit <- stan(
  file    = "classroom_ri.stan",
  data    = stan_data,
  chains  = 2,
  iter    = 6000,
  warmup  = 1000,
  seed    = 123,
  control = list(adapt_delta = 0.95)
)

print(fit, pars = c("beta", "sigma_b", "sigma_e", "icc"))
