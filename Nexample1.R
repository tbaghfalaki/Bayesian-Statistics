rm(list=ls())
library(sn)
library(R2jags)

data(ais)

# Extract variables
bmi <- log(ais$BMI)
sex <- ifelse(ais$sex == "female", 0, 1)

# Data for JAGS
data_jags <- list(
  N = length(bmi),
  bmi = bmi,
  sex = sex
)

# Initial values
inits_fun <- function() {
  list(
    beta0 = mean(bmi),
    beta1 = 0,
    tau2 = 1 / var(bmi)
  )
}

# Parameters to monitor
params <- c("beta0", "beta1", "sigma2")

# JAGS model
linear_model <- "
model {

  for (i in 1:N) {
    bmi[i] ~ dnorm(mu[i], tau2)
    mu[i] <- beta0 + beta1 * sex[i]
  }

  # Priors
  beta0 ~ dnorm(0, 0.0001)
  beta1 ~ dnorm(0, 0.0001)

  # Precision
  tau2 ~ dgamma(0.001, 0.001)

  # Standard deviation
  sigma2 <- 1 / sqrt(tau2)
}
"

# Fit Bayesian regression model
fit <- jags(
  data = data_jags,
  inits = inits_fun,
  parameters.to.save = params,
  model.file = textConnection(linear_model),
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 2000
)

# Results
print(fit)

# Posterior summaries
fit$summary

# Trace plots
plot(fit)