# Load libraries
library(rjags)
library(coda)
library(sn)
data(ais)

# -----------------------
# Prepare data
# -----------------------
ais$sex_num <- ifelse(ais$sex == "male", 1, 0) # 0=female, 1=male
ais$Ht_s <- (ais$Ht - mean(ais$Ht))/sd(ais$Ht)
ais$Wt_s <- (ais$Wt - mean(ais$Wt))/sd(ais$Wt)
y <- ais$BMI
N <- length(y)

# -----------------------
# 1. Model 1: BMI ~ sex
# -----------------------
x1 <- ais$sex_num
data1 <- list(y=y, x=x1, N=N)

model1_string <- "
model {
  for(i in 1:N){
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

inits1 <- function() list(beta0=0, beta1=0, tau=1)
parameters1 <- c("beta0","beta1","sigma")

mod1 <- jags.model("model1.txt", data=data1, inits=inits1, n.chains=3, quiet=TRUE)
update(mod1, 1000) # burn-in
samp1 <- coda.samples(mod1, variable.names=parameters1, n.iter=5000)
summary(samp1)
# Posterior samples
samp1_mat <- as.matrix(samp1)
beta0_s <- samp1_mat[,"beta0"]
beta1_s <- samp1_mat[,"beta1"]
sigma_s <- samp1_mat[,"sigma"]

# Compute DIC
dic1 <- dic.samples(mod1, n.iter=5000, type="pD")
dic1_val <- sum(dic1$deviance) + sum(dic1$penalty)

# LPML
cpo <- numeric(N)
S <- length(beta0_s)
for(i in 1:N){
  lik <- dnorm(y[i], mean = beta0_s + beta1_s * x1[i], sd = sigma_s)
  cpo[i] <- 1 / mean(1/lik)
}
lpml1 <- sum(log(cpo))

# Approximate BIC
beta0_hat <- mean(beta0_s)
beta1_hat <- mean(beta1_s)
sigma_hat <- mean(sigma_s)
loglik <- sum(dnorm(y, mean = beta0_hat + beta1_hat*x1, sd = sigma_hat, log=TRUE))
bic1 <- 3*log(N) - 2*loglik  # 2 betas + sigma

# -----------------------
# 2. Model 2: BMI ~ sex + height + weight
# -----------------------
x2_sex <- ais$sex_num
x2_Ht <- ais$Ht_s
x2_Wt <- ais$Wt_s
data2 <- list(y=y, x1=x2_sex, x2=x2_Ht, x3=x2_Wt, N=N)

model2_string <- "
model {
  for(i in 1:N){
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

inits2 <- function() list(beta0=0, beta1=0, beta2=0, beta3=0, tau=1)
parameters2 <- c("beta0","beta1","beta2","beta3","sigma")

mod2 <- jags.model("model2.txt", data=data2, inits=inits2, n.chains=3, quiet=TRUE)
update(mod2, 1000) # burn-in
samp2 <- coda.samples(mod2, variable.names=parameters2, n.iter=5000)
summary(samp2)
# Posterior samples
samp2_mat <- as.matrix(samp2)
beta0_s <- samp2_mat[,"beta0"]
beta1_s <- samp2_mat[,"beta1"]
beta2_s <- samp2_mat[,"beta2"]
beta3_s <- samp2_mat[,"beta3"]
sigma_s <- samp2_mat[,"sigma"]

# Compute DIC
dic2 <- dic.samples(mod2, n.iter=5000, type="pD")
dic2_val <- sum(dic2$deviance) + sum(dic2$penalty)

# LPML
cpo <- numeric(N)
S <- length(beta0_s)
for(i in 1:N){
  lik <- dnorm(y[i], mean = beta0_s + beta1_s*x2_sex[i] + beta2_s*x2_Ht[i] + beta3_s*x2_Wt[i], sd = sigma_s)
  cpo[i] <- 1 / mean(1/lik)
}
lpml2 <- sum(log(cpo))

# Approximate BIC
beta0_hat <- mean(beta0_s)
beta1_hat <- mean(beta1_s)
beta2_hat <- mean(beta2_s)
beta3_hat <- mean(beta3_s)
sigma_hat <- mean(sigma_s)
loglik <- sum(dnorm(y, mean = beta0_hat + beta1_hat*x2_sex + beta2_hat*x2_Ht + beta3_hat*x2_Wt, sd = sigma_hat, log=TRUE))
bic2 <- 5*log(N) - 2*loglik  # 4 betas + sigma


