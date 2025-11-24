rm(list=ls())

# Libraries
library(sn)       # AIS dataset
library(MASS)     # mvrnorm
library(ggplot2)
library(dplyr)

# Prepare data
data(ais)
ais$sex_num <- ifelse(ais$sex == "male", 1, 0) # 0=female, 1=male

y <- as.matrix(ais$BMI)          # N x 1
N <- nrow(y)
X <- as.matrix(cbind(1, ais$sex_num))  # intercept + sex
d <- ncol(X)

# Prior parameters
tau2 <- 100
a0 <- 0.1
b0 <- 0.1

# Initialize variational parameters
mu <- matrix(0, nrow=d, ncol=1)         # d x 1
Lambda <- diag(1/tau2, d)               # d x d (precision)
a_tilde <- a0 + N/2
b_tilde <- b0 + 0.5 * sum((y - X %*% mu)^2)

# CAVI updates
max_iter <- 1000
tol <- 1e-6

for(iter in 1:max_iter){
  
  # Ensure scalar
  scal <- as.numeric(a_tilde / b_tilde)
  
  # Update q(beta)
  Lambda_new <- scal * t(X) %*% X + diag(1/tau2, d)
  mu_new <- solve(Lambda_new) %*% (scal * t(X) %*% y)
  
  # Update q(sigma^2)
  quad <- t(y - X %*% mu_new) %*% (y - X %*% mu_new)    # 1x1 matrix
  trace_term <- sum(diag(X %*% solve(Lambda_new) %*% t(X)))
  b_tilde_new <- b0 + 0.5 * (as.numeric(quad) + trace_term)
  
  # Convergence check
  if(max(abs(mu - mu_new)) < tol && max(abs(diag(Lambda - Lambda_new))) < tol && abs(b_tilde - b_tilde_new) < tol){
    mu <- mu_new
    Lambda <- Lambda_new
    b_tilde <- b_tilde_new
    break
  }
  
  # Update
  mu <- mu_new
  Lambda <- Lambda_new
  b_tilde <- b_tilde_new
}

# Posterior covariance of beta
Sigma_beta <- solve(Lambda)   # covariance matrix
sd_beta <- sqrt(diag(Sigma_beta))

# Posterior mean and variance for sigma^2 (if defined)
mean_sigma2 <- if(a_tilde > 1) b_tilde / (a_tilde - 1) else NA
var_sigma2  <- if(a_tilde > 2) (b_tilde^2) / ((a_tilde - 1)^2 * (a_tilde - 2)) else NA

# Print results
cat("Variational posterior of beta (Normal):\n")
cat("Mean (mu):\n"); print(as.vector(mu))
cat("\nCovariance (Sigma_beta):\n"); print(Sigma_beta)
cat("\nStandard deviations (sd_beta):\n"); print(sd_beta)

cat("\nVariational posterior of sigma^2 (Inverse-Gamma):\n")
cat("a_tilde =", a_tilde, "\n")
cat("b_tilde =", b_tilde, "\n")
cat("E[sigma^2] =", mean_sigma2, "\n")
cat("Var[sigma^2] =", var_sigma2, "\n")



# Posterior samples of betas and sigma2
nsamp <- 10000
beta_samples <- mvrnorm(nsamp, mu = as.vector(mu), Sigma = Sigma_beta)
# sample sigma^2 from inverse-gamma via 1 / Gamma(shape=a_tilde, rate=b_tilde)
sigma2_samples <- 1 / rgamma(nsamp, shape = a_tilde, rate = b_tilde)

# Prepare long dataframe for ggplot (beta0, beta1, sigma2)
df <- data.frame(
  value = c(beta_samples[,1], beta_samples[,2], sigma2_samples),
  parameter = rep(c("beta0", "beta1", "sigma2"), each = nsamp)
)

# Faceted density plot with separate fills
p_facet <- ggplot(df, aes(x = value, fill = parameter)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~parameter, scales = "free") +
  scale_fill_manual(values = c("skyblue", "lightgreen", "orange")) +
  labs(title = "VB Posterior Densities (AIS: BMI ~ sex)",
       x = "Value", y = "Density") +
  theme_bw() +
  theme(legend.position = "none")

print(p_facet)

 
