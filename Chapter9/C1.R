rm(list=ls())
library(R2jags)

# Load data
Data <- read.csv("C:\\Users\\p80744tb\\Desktop\\titanic\\train_complete.csv", stringsAsFactors = FALSE)
names(Data)
attach(Data)
head(Data)

 

# Create design matrix (includes intercept)
X <- model.matrix(~ Sex + Age + as.factor(Pclass) + Fare + 
                    as.factor(Embarked) + SibSp + Parch, data = Data)
head(X)
 

K <- ncol(X)


# Prepare data list for JAGS
data_jags <- list(
  N = nrow(X),
  K = K,
  Survived = Data$Survived,
  X = X
)


# Initial values function
inits_fun <- function() {
  list(beta = rep(0, K))   # initialize all betas to 0
}


# Parameters to monitor
params <- c("beta")


# JAGS model as a text string
logistic_model <- "
model {
  # Likelihood
  for (i in 1:N) {
    Survived[i] ~ dbern(pi[i])
    logit(pi[i]) <- inprod(beta[], X[i, ])
  }

  # Priors
  for (k in 1:K) {
    beta[k] ~ dnorm(0, 0.001)
  }
}
"


# Fit the model using R2jags
fit <- jags(
  data = data_jags,
  inits = inits_fun,
  parameters.to.save = params,
  model.file = textConnection(logistic_model),
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 2000
)


# Summarize results
print(fit)




