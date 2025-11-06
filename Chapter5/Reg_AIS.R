# Load required packages
library(sn)

# Load AIS data
data(ais)
attach(ais)

# Encode sex as binary: 0 = female, 1 = male
ais$sex_num <- ifelse(ais$sex == "male", 1, 0)

set.seed(123)
n <- nrow(ais)
y <- ais$BMI
x <- ais$sex_num
n_train <- length(y)

# Set hyperparameters
mu_beta <- 0
sigma_beta2 <- 1
a_sigma <- 0.01
b_sigma <- 0.01

# Gibbs sampler settings
n_iter <- 10000
burn_in <- 5000

beta0_samples <- numeric(n_iter)
beta1_samples <- numeric(n_iter)
sigma2_samples <- numeric(n_iter)

# Initialize
beta0 <- 0
beta1 <- 0
sigma2 <- 1

# Gibbs sampler
for (iter in 1:n_iter) {
  
  # Update beta0 conditional on beta1, sigma2, y
  sigma_beta0_post <- 1 / (n_train / sigma2 + 1 / sigma_beta2)
  mu_beta0_post <- sigma_beta0_post * (sum(y - beta1 * x) / sigma2 + mu_beta / sigma_beta2)
  beta0 <- rnorm(1, mean = mu_beta0_post, sd = sqrt(sigma_beta0_post))
  
  # Update beta1 conditional on beta0, sigma2, y
  sigma_beta1_post <- 1 / (sum(x^2) / sigma2 + 1 / sigma_beta2)
  mu_beta1_post <- sigma_beta1_post * (sum(x * (y - beta0)) / sigma2 + mu_beta / sigma_beta2)
  beta1 <- rnorm(1, mean = mu_beta1_post, sd = sqrt(sigma_beta1_post))
  
  # Update sigma2 conditional on beta0, beta1, y
  a_post <- a_sigma + n_train / 2
  b_post <- b_sigma + sum((y - beta0 - beta1 * x)^2) / 2
  sigma2 <- 1 / rgamma(1, shape = a_post, rate = b_post)  # Inverse-Gamma
  
  # Store samples
  beta0_samples[iter] <- beta0
  beta1_samples[iter] <- beta1
  sigma2_samples[iter] <- sigma2
}

# Apply burn-in
beta0_post <- beta0_samples[(burn_in + 1):n_iter]
beta1_post <- beta1_samples[(burn_in + 1):n_iter]
sigma2_post <- sigma2_samples[(burn_in + 1):n_iter]

# Posterior summaries
posterior_summary <- data.frame(
  Parameter = c("beta0", "beta1", "sigma2"),
  Mean = c(mean(beta0_post), mean(beta1_post), mean(sigma2_post)),
  SD = c(sd(beta0_post), sd(beta1_post), sd(sigma2_post)),
  `2.5%` = c(quantile(beta0_post, 0.025), quantile(beta1_post, 0.025), quantile(sigma2_post, 0.025)),
  `97.5%` = c(quantile(beta0_post, 0.975), quantile(beta1_post, 0.975), quantile(sigma2_post, 0.975))
)

print(posterior_summary)

# Trace plots with Greek letters
par(mfrow = c(3,1))
plot(beta0_post, type='l', main=expression(paste("Trace plot: ", beta[0])), ylab=expression(beta[0]))
plot(beta1_post, type='l', main=expression(paste("Trace plot: ", beta[1])), ylab=expression(beta[1]))
plot(sigma2_post, type='l', main=expression(paste("Trace plot: ", sigma^2)), ylab=expression(sigma^2))

# Posterior density plots
par(mfrow = c(3,1))
plot(density(beta0_post), main=expression(paste("Posterior density: ", beta[0])), xlab=expression(beta[0]))
plot(density(beta1_post), main=expression(paste("Posterior density: ", beta[1])), xlab=expression(beta[1]))
plot(density(sigma2_post), main=expression(paste("Posterior density: ", sigma^2)), xlab=expression(sigma^2))




# Improved trace and posterior density plots
par(mfrow = c(3,2), mar = c(4,4,2,1))  # 3 rows, 2 columns, adjust margins

# Beta0
plot(beta0_post, type='l', col='orange3', lwd=1.5,
     main=expression(paste("Trace plot: ", beta[0])),
     ylab=expression(beta[0]), xlab="Iteration")
plot(density(beta0_post), col='red', lwd=2,
     main=expression(paste("Posterior density: ", beta[0])),
     xlab=expression(beta[0]))
polygon(density(beta0_post), col=rgb(1,0,0,0.4), border=NA)

# Beta1
plot(beta1_post, type='l', col='orange3', lwd=1.5,
     main=expression(paste("Trace plot: ", beta[1])),
     ylab=expression(beta[1]), xlab="Iteration")
plot(density(beta1_post), col='red', lwd=2,
     main=expression(paste("Posterior density: ", beta[1])),
     xlab=expression(beta[1]))
polygon(density(beta1_post), col=rgb(1,0,0,0.4), border=NA)

# Sigma^2
plot(sigma2_post, type='l', col='orange3', lwd=1.5,
     main=expression(paste("Trace plot: ", sigma^2)),
     ylab=expression(sigma^2), xlab="Iteration")
plot(density(sigma2_post), col='red', lwd=2,
     main=expression(paste("Posterior density: ", sigma^2)),
     xlab=expression(sigma^2))
polygon(density(sigma2_post), col=rgb(1,0,0,0.4), border=NA)
