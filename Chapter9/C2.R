# Load libraries
library(R2jags)
library(bellreg)
library(dplyr)

# Load data
data(cells)
attach(cells)

# Prepare dataset
cells <- cells %>%
  mutate(
    smoker = as.numeric(smoker),
    gender = as.numeric(gender),
    weightf = factor(weight, levels = c("normal", "over", "obese")),
    agef = factor(age, levels = c("young", "mid", "old"))
  )

# Model matrix: Intercept, Smoking, BMI(Overweight/Obese), Age(Middle/Old), Gender
X <- model.matrix(~ smoker + weightf + agef + gender, data = cells)
table(X[,3])
K <- ncol(X)
N <- nrow(X)
Y <- cells$cells   # response variable (count of infected cells)

# Prepare data list for JAGS
data_jags <- list(
  N = N,
  K = K,
  Y = Y,
  X = X
)

# Initial values
inits_fun <- function() {
  list(beta = rep(0, K))
}

# Parameters to monitor
params <- c("beta")

# JAGS model (Poisson)
poisson_model <- "
model {
  # Likelihood
  for (i in 1:N) {
    Y[i] ~ dpois(mu[i])
    log(mu[i]) <- inprod(beta[], X[i, ])
  }

  # Priors
  for (k in 1:K) {
    beta[k] ~ dnorm(0, 0.001)
  }
}
"

# Fit Poisson model
fit_pois <- jags(
  data = data_jags,
  inits = inits_fun,
  parameters.to.save = params,
  model.file = textConnection(poisson_model),
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 2000
)

# Summarize results
print(fit_pois)



# JAGS model (Negative Binomial)
nb_model <- "
model {
  for (i in 1:N) {
    # Poisson-Gamma mixture representation
    Y[i] ~ dpois(lambda[i])
    lambda[i] ~ dgamma(r, r / mu[i])
    log(mu[i]) <- inprod(beta[], X[i, ])
  }

  # Priors for regression coefficients
  for (k in 1:K) {
    beta[k] ~ dnorm(0, 0.001)
  }

  # Prior for overdispersion parameter r
  r ~ dunif(0, 100)
}
"

# Updated parameters to monitor (include r)
params_nb <- c("beta", "r")

# Initial values
inits_nb <- function() {
  list(beta = rep(0, K), r = runif(1, 1, 20))
}

# Fit Negative Binomial model
fit_nb <- jags(
  data = data_jags,
  inits = inits_nb,
  parameters.to.save = params_nb,
  model.file = textConnection(nb_model),
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 2000
)

# Summarize results
print(fit_pois)
print(fit_nb)



fit_pois$BUGSoutput$DIC
fit_nb$BUGSoutput$DIC
