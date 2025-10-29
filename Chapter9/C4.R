# Load libraries
library(R2jags)


# Simulate data
set.seed(123)
n <- 100
y <- rnorm(n, mean = 0, sd = 1)
zeros <- rep(0, n)   # dummy variable for zeros-trick
C <- 10              # large constant for zeros-trick


# Standard JAGS model
data_standard <- list(y = y, n = n)

model_standard <- "
model {
  for (i in 1:n) {
    y[i] ~ dnorm(mu, tau)
  }
  mu ~ dnorm(0, 0.001)
  tau <- 1/pow(sigma,2)
  sigma ~ dunif(0, 10)
}
"

fit_standard <- jags(
  data = data_standard,
  parameters.to.save = c("mu", "sigma"),
  model.file = textConnection(model_standard),
  n.chains = 3,
  n.iter = 5000,
  n.burnin = 1000
)

print("Standard JAGS Results:")
print(fit_standard)


# Zeros-Trick JAGS model
data_zeros <- list(y = y, zeros = zeros, n = n, C = C)

model_zeros <- "
model {
  for (i in 1:n) {
    zeros[i] ~ dpois(lambda[i])
    lambda[i] <- -log(1/sqrt(2*3.14159*pow(sigma,2)) * 
    exp(-pow((y[i]-mu)/sigma,2)/2)) + C
  }

  # Priors
  mu ~ dnorm(0, 0.001)
  sigma ~ dunif(0, 100)
}
"

fit_zeros <- jags(
  data = data_zeros,
  parameters.to.save = c("mu", "sigma"),
  model.file = textConnection(model_zeros),
  n.chains = 3,
  n.iter = 5000,
  n.burnin = 1000
)

print("Zeros-Trick JAGS Results:")
print(fit_zeros)
