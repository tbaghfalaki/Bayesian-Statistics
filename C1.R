# Load required libraries
library(MASS)      # for mvrnorm
library(ggplot2)   # for plotting
library(reshape2)  # for melting matrices

# AR(1) parameters
n <- 20
rho <- 0.8
sigma2 <- 1

# Construct covariance matrix for AR(1)
Sigma <- matrix(0, n, n)
for(i in 1:n){
  for(j in 1:n){
    Sigma[i, j] <- sigma2 * rho^abs(i-j) / (1 - rho^2)
  }
}

# Compute precision matrix (inverse covariance)
Q <- solve(Sigma)

# Simulate one realization of AR(1)
set.seed(123)
eta <- mvrnorm(1, mu = rep(0, n), Sigma = Sigma)

# Plot the latent field
plot(1:n, eta, type="o", pch=16, col="blue",
     xlab="Index i", ylab=expression(eta[i]),
     main="AR(1) Process as a GMRF")
grid()

# Plot ACF of the latent field
acf(eta, main="ACF of AR(1) Latent Field", lag.max=20,lwd=2)

# Heatmap of covariance matrix
Sigma_melt <- melt(Sigma)
colnames(Sigma_melt) <- c("i", "j", "Covariance")

ggplot(Sigma_melt, aes(x=i, y=j, fill=Covariance)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0) +
  labs(title="Covariance Matrix Heatmap", x="Index i", y="Index j") +
  theme_minimal()

# Heatmap of precision matrix
Q_melt <- melt(Q)
colnames(Q_melt) <- c("i", "j", "Precision")

ggplot(Q_melt, aes(x=i, y=j, fill=Precision)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0) +
  labs(title="Precision Matrix Heatmap (Sparsity Pattern)", x="Index i", y="Index j") +
  theme_minimal()
