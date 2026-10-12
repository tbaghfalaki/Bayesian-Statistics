rm(list=ls())
# Load libraries
library(R2jags)
library(coda)
library(sn)     # for the AIS dataset


# Prepare data
data(ais)
ais$sex_num <- ifelse(ais$sex == "male", 1, 0) # 0=female, 1=male
ais$Ht_s <- (ais$Ht - mean(ais$Ht))/sd(ais$Ht)
ais$Wt_s <- (ais$Wt - mean(ais$Wt))/sd(ais$Wt)
y <- ais$BMI
N <- length(y)


# Model 1: BMI ~ sex
model1_string <- "
model {
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 * x[i]
  }
  beta0 ~ dnorm(0.0, 1.0E-6)
  beta1 ~ dnorm(0.0, 1.0E-6)
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)
}
"
writeLines(model1_string, "model1.txt")

# Data and parameters
x1 <- ais$sex_num
data1 <- list(y = y, x = x1, N = N)
inits1 <- function() list(beta0 = 0, beta1 = 0, tau = 1)
params1 <- c("beta0", "beta1", "sigma")


# Run Model 1 with R2jags

mod1 <- jags(
  data = data1,
  inits = inits1,
  parameters.to.save = params1,
  model.file = "model1.txt",
  n.chains = 3,
  n.iter = 6000,
  n.burnin = 1000,
  DIC = TRUE
)

print(mod1)
samp1 <- as.mcmc(mod1)
summary(samp1)

# Posterior samples
samp1_mat <- mod1$BUGSoutput$sims.matrix
beta0_s <- samp1_mat[,"beta0"]
beta1_s <- samp1_mat[,"beta1"]
sigma_s <- samp1_mat[,"sigma"]


# Compute model metrics


# DIC
dic1_val <- mod1$BUGSoutput$DIC

# LPML (Conditional Predictive Ordinates)
cpo <- numeric(N)
S <- length(beta0_s)
for (i in 1:N) {
  lik <- dnorm(y[i], mean = beta0_s + beta1_s * x1[i], sd = sigma_s)
  cpo[i] <- 1 / mean(1 / lik)
}
lpml1 <- sum(log(cpo))

# Approximate BIC
beta0_hat <- mean(beta0_s)
beta1_hat <- mean(beta1_s)
sigma_hat <- mean(sigma_s)
loglik <- sum(dnorm(y, mean = beta0_hat + beta1_hat * x1, sd = sigma_hat, log = TRUE))
bic1 <- 3 * log(N) - 2 * loglik  # 2 betas + sigma



# Model 2: BMI ~ sex + height + weight
model2_string <- "
model {
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- beta0 + beta1*x1[i] + beta2*x2[i] + beta3*x3[i]
  }
  beta0 ~ dnorm(0.0, 1.0E-6)
  beta1 ~ dnorm(0.0, 1.0E-6)
  beta2 ~ dnorm(0.0, 1.0E-6)
  beta3 ~ dnorm(0.0, 1.0E-6)
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)
}
"
writeLines(model2_string, "model2.txt")

# Data and parameters
x2_sex <- ais$sex_num
x2_Ht <- ais$Ht_s
x2_Wt <- ais$Wt_s
data2 <- list(y = y, x1 = x2_sex, x2 = x2_Ht, x3 = x2_Wt, N = N)
inits2 <- function() list(beta0 = 0, beta1 = 0, beta2 = 0, beta3 = 0, tau = 1)
params2 <- c("beta0", "beta1", "beta2", "beta3", "sigma")


# Run Model 2 with R2jags
mod2 <- jags(
  data = data2,
  inits = inits2,
  parameters.to.save = params2,
  model.file = "model2.txt",
  n.chains = 3,
  n.iter = 6000,
  n.burnin = 1000,
  DIC = TRUE
)

print(mod2)
samp2 <- as.mcmc(mod2)
summary(samp2)

# Posterior samples
samp2_mat <- mod2$BUGSoutput$sims.matrix
beta0_s <- samp2_mat[,"beta0"]
beta1_s <- samp2_mat[,"beta1"]
beta2_s <- samp2_mat[,"beta2"]
beta3_s <- samp2_mat[,"beta3"]
sigma_s <- samp2_mat[,"sigma"]


# Compute model metrics


# DIC
dic2_val <- mod2$BUGSoutput$DIC

# LPML
cpo <- numeric(N)
S <- length(beta0_s)
for (i in 1:N) {
  lik <- dnorm(y[i], mean = beta0_s + beta1_s*x2_sex[i] + beta2_s*x2_Ht[i] + beta3_s*x2_Wt[i], sd = sigma_s)
  cpo[i] <- 1 / mean(1 / lik)
}
lpml2 <- sum(log(cpo))

# Approximate BIC
beta0_hat <- mean(beta0_s)
beta1_hat <- mean(beta1_s)
beta2_hat <- mean(beta2_s)
beta3_hat <- mean(beta3_s)
sigma_hat <- mean(sigma_s)
loglik <- sum(dnorm(y, mean = beta0_hat + beta1_hat*x2_sex + beta2_hat*x2_Ht + beta3_hat*x2_Wt, sd = sigma_hat, log=TRUE))
bic2 <- 5 * log(N) - 2 * loglik  # 4 betas + sigma



# Model 1
cat("DIC:", dic1_val, "\nLPML:", lpml1, "\nBIC:", bic1, "\n")
# Model 2
cat("DIC:", dic2_val, "\nLPML:", lpml2, "\nBIC:", bic2, "\n")
