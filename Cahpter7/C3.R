rm(list=ls())

# --- Libraries ---
library(sn)
library(MASS)
library(ggplot2)
library(dplyr)

# --- Prepare data ---
data(ais)
ais$sex_num <- ifelse(ais$sex == "male", 1, 0)

y <- as.matrix(ais$BMI)
N <- nrow(y)
X <- as.matrix(cbind(1, ais$sex_num))
d <- ncol(X)

# --- Prior parameters ---
tau2 <- 100
a0 <- 0.01
b0 <- 0.01

# --- Initialize VB parameters ---
mu <- matrix(0, nrow=d, ncol=1)
Lambda <- diag(1/tau2, d)
a_tilde <- a0 + N/2
b_tilde <- b0 + 0.5 * sum((y - X %*% mu)^2)

# --- CAVI updates ---
max_iter <- 1000
tol <- 0.01

for(iter in 1:max_iter){
  scal <- as.numeric(a_tilde / b_tilde)
  Lambda_new <- scal * t(X) %*% X + diag(1/tau2, d)
  mu_new <- solve(Lambda_new) %*% (scal * t(X) %*% y)
  
  quad <- t(y - X %*% mu_new) %*% (y - X %*% mu_new)
  trace_term <- sum(diag(X %*% solve(Lambda_new) %*% t(X)))
  b_tilde_new <- b0 + 0.5 * (as.numeric(quad) + trace_term)
  
  if(max(abs(mu - mu_new)) < tol && max(abs(diag(Lambda - Lambda_new))) < tol && abs(b_tilde - b_tilde_new) < tol){
    mu <- mu_new
    Lambda <- Lambda_new
    b_tilde <- b_tilde_new
    break
  }
  
  mu <- mu_new
  Lambda <- Lambda_new
  b_tilde <- b_tilde_new
}

# --- Posterior covariance of beta ---
Sigma_beta <- solve(Lambda)

# --- Posterior samples ---
nsamp <- 1000
beta_samples <- mvrnorm(nsamp, mu = as.vector(mu), Sigma = Sigma_beta)
sigma2_samples <- 1 / rgamma(nsamp, shape = a_tilde, rate = b_tilde)

# --- Exact posterior for overlay ---
y_bar <- mean(y)
XtX <- t(X) %*% X
beta_hat <- solve(XtX) %*% t(X) %*% y

a_n <- a0 + N/2
b_n <- b0 + 0.5 * sum((y - X %*% beta_hat)^2) + (t(beta_hat) %*% beta_hat) / (2*tau2)

theta_seq <- seq(min(beta_samples[,1])-1, max(beta_samples[,1])+1, length.out=500)
theta2_seq <- seq(min(beta_samples[,2])-1, max(beta_samples[,2])+1, length.out=500)
sigma2_seq <- seq(min(sigma2_samples)*0.5, max(sigma2_samples)*1.5, length.out=500)

# Corrected inverse-gamma density function
dinvgamma <- function(x, shape, rate){
  x <- as.vector(x)
  (rate^shape / gamma(shape)) * x^(-shape-1) * exp(-rate/x)
}

# Prepare exact density dataframes
df_exact <- rbind(
  data.frame(value=theta_seq, density=dnorm(theta_seq, mean=beta_hat[1], sd=sqrt(b_n/(a_n*XtX[1,1]))),
             parameter="beta0", type="Exact"),
  data.frame(value=theta2_seq, density=dnorm(theta2_seq, mean=beta_hat[2], sd=sqrt(b_n/(a_n*XtX[2,2]))),
             parameter="beta1", type="Exact"),
  data.frame(value=sigma2_seq, density=dinvgamma(sigma2_seq, a_n, b_n),
             parameter="sigma2", type="Exact")
)

# Prepare VB shaded dataframes
df_vb <- rbind(
  data.frame(value=beta_samples[,1], parameter="beta0", type="VB"),
  data.frame(value=beta_samples[,2], parameter="beta1", type="VB"),
  data.frame(value=sigma2_samples, parameter="sigma2", type="VB")
)

# --- Plot ---
# Optional: add dashed lines for estimated beta0, beta1, sigma2
true_vals <- data.frame(
  parameter=c("beta0", "beta1", "sigma2"),
  true=c(beta_hat[1], beta_hat[2], mean(sigma2_samples))  # or true values if known
)

ggplot() +
  geom_histogram(data=subset(df_vb, type=="VB"), aes(x=value, fill=parameter, y=..density..),
                 bins=60, alpha=0.4, position="identity") +
  geom_line(data=df_exact, aes(x=value, y=density, color=parameter), size=0.5) +
  geom_vline(data=true_vals, aes(xintercept=true, color=parameter), linetype="dashed", size=0.5) +
  facet_wrap(~parameter, scales="free") +
  scale_fill_manual(values=c("beta0"="lightgreen","beta1"="skyblue","sigma2"="orange")) +
  scale_color_manual(values=c("beta0"="darkgreen","beta1"="blue","sigma2"="red")) +
  labs(title="",
       x="Value", y="Density") +
  theme_bw() +
  theme(legend.position="none")
