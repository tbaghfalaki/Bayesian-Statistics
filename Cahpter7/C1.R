### ---- Variational Bayes for Normal Model (Example 2) ----
rm(list=ls())
set.seed(9)

# --- Generate data ---
n <- 100
theta_true <- 2
sigma2_true <- 1
y <- rnorm(n, mean = theta_true, sd = sqrt(sigma2_true))

# --- Prior parameters ---
tau2 <- 100         # prior variance for theta
a0 <- 0.01          # IG prior shape
b0 <- 0.01          # IG prior scale

# --- Initialize VB parameters ---
mu <- 0
lambda <- 1
tilde_a <- a0
tilde_b <- b0

max_iter <- 2000
tol <- 1e-8

# --- Coordinate Ascent Variational Inference ---
for (iter in 1:max_iter) {
  
  mu_old <- mu
  lambda_old <- lambda
  tb_old <- tilde_b
  
  # update q(theta)
  lambda <- n * (tilde_a / tilde_b) + 1/tau2
  mu <- ( (tilde_a / tilde_b) * sum(y) ) / lambda
  
  # update q(sigma^2)
  tilde_a <- a0 + n/2
  tilde_b <- b0 + 0.5 * sum( (y - mu)^2 + 1/lambda )
  
  # convergence check
  if (max(abs(mu - mu_old), abs(lambda - lambda_old), abs(tilde_b - tb_old)) < tol)
    break
}
print(iter)

# Posterior mean and SD analytically
theta_mean_vb <- mu
theta_sd_vb   <- sqrt(1/lambda)

sigma2_mean_vb <- tilde_b / (tilde_a - 1)
sigma2_sd_vb   <- sqrt(tilde_b^2 / ((tilde_a - 1)^2 * (tilde_a - 2)))

cat("VB posterior mean of theta:", theta_mean_vb, "\n")
cat("VB posterior SD of theta:", theta_sd_vb, "\n\n")
cat("VB posterior mean of sigma^2:", sigma2_mean_vb, "\n")
cat("VB posterior SD of sigma^2:", sigma2_sd_vb, "\n")

# --- Draw posterior samples from VB ---
S <- 1000
theta_samples <- rnorm(S, mean = mu, sd = sqrt(1/lambda))
sigma2_samples <- 1 / rgamma(S, shape = tilde_a, rate = tilde_b)

# --- Exact posterior densities ---
theta_seq <- seq(1.5, 2.5, length.out=500)
sigma2_seq <- seq( 0.5, max(sigma2_samples)*1.5, length.out=500)

# Inverse-Gamma density
dinvgamma <- function(x, shape, rate){
  (rate^shape / gamma(shape)) * x^(-shape-1) * exp(-rate/x)
}

# Exact posterior parameters
a_n <- a0 + n/2
b_n <- b0 + 0.5 * sum((y - mean(y))^2) + (n*(mean(y) - 0)^2)/(2*(tau2+n))

# Exact densities
sigma2_density <- dinvgamma(sigma2_seq, a_n, b_n)
theta_df <- 2 * a_n
theta_mean <- mean(y)
theta_scale <- sqrt(b_n / (a_n * n))
theta_density <- dt((theta_seq - theta_mean)/theta_scale, df=theta_df)/theta_scale

# --- Prepare data frames for plotting ---
df_vb <- data.frame(
  value = c(theta_samples, sigma2_samples),
  parameter = rep(c("theta","sigma2"), each=S),
  type = "VB"
)

df_exact <- data.frame(
  value = c(theta_seq, sigma2_seq),
  density = c(theta_density, sigma2_density),
  parameter = rep(c("theta","sigma2"), each=length(theta_seq)),
  type = "Exact"
)

# True values
true_vals <- data.frame(
  parameter = c("theta","sigma2"),
  true = c(theta_true, sigma2_true)
)

# --- Plot VB (shaded) vs Exact (line) ---
ggplot() +
  geom_density(data=df_vb, aes(x=value, fill=type), alpha=0.5, color=NA) +
  geom_line(data=df_exact, aes(x=value, y=density, color=type), size=1.2) +
  geom_vline(data=true_vals, aes(xintercept=true), linetype="dashed", size=1) +
  facet_wrap(~parameter, scales="free") +
  scale_fill_manual(values=c("VB"="skyblue","Exact"="pink")) +
  scale_color_manual(values=c("VB"="skyblue","Exact"="red")) +
  labs(title="VB Posterior (shaded) vs Exact Posterior (line)",
       x="Value", y="Density", fill="Distribution", color="Distribution") +
  theme_bw()
