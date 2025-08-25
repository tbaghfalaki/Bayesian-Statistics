# Load required packages
library(sn)
library(rjags)
library(coda)
library(ggplot2)

# Load AIS data
data(ais)
attach(ais)

# Show variable names
names(ais)

# Encode sex as binary: 0 = female, 1 = male
ais$sex_num <- ifelse(ais$sex == "male", 1, 0)

# Split data: 70% training, 30% validation
set.seed(123)
n <- nrow(ais)
train_indices <- sample(1:n, size = round(0.7 * n))
train_data <- ais[train_indices, ]
valid_data <- ais[-train_indices, ]

# Prepare data for JAGS
jags_data <- list(
  y = train_data$BMI,
  x = train_data$sex_num,
  N = nrow(train_data)
)

# --------- JAGS Model 1: Normal Errors ---------
model_normal <- "
model {
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + beta * x[i]
  }
  alpha ~ dnorm(0, 0.001)
  beta ~ dnorm(0, 0.001)
  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 100)
}
"

# Fit normal model
jags_normal <- jags.model(textConnection(model_normal), data = jags_data, n.chains = 3)
update(jags_normal, 1000)  # burn-in
samples_normal <- coda.samples(jags_normal, c("alpha", "beta", "sigma"), 5000)

# --------- JAGS Model 2: Student-t Errors ---------
model_t <- "
model {
  for (i in 1:N) {
    y[i] ~ dt(mu[i], tau, nu)
    mu[i] <- alpha + beta * x[i]
  }
  alpha ~ dnorm(0, 0.001)
  beta ~ dnorm(0, 0.001)
  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 100)
  nu ~ dunif(1, 30)
}
"

# Fit t model
jags_t <- jags.model(textConnection(model_t), data = jags_data, n.chains = 3)
update(jags_t, 1000)  # burn-in
samples_t <- coda.samples(jags_t, c("alpha", "beta", "sigma", "nu"), 5000)

# --------- Extract Posterior Summaries ---------
mean_normal <- summary(samples_normal)$statistics[, "Mean"]
mean_t <- summary(samples_t)$statistics[, "Mean"]



summary_normal <- summary(samples_normal)$statistics 
summary_t <- summary(samples_t)$statistics 
#xtable::xtable(cbind(summary_normal[,1:2]))
#xtable::xtable(cbind(summary_t[,1:2]))

# --------- Histogram and Overlay of Densities (Males) ---------

library(ggplot2)
library(patchwork)  # install.packages("patchwork") if needed

# Predicted means
mu_female_normal <- mean_normal["alpha"]
mu_female_t <- mean_t["alpha"]

mu_male_normal <- mean_normal["alpha"] + mean_normal["beta"]
mu_male_t <- mean_t["alpha"] + mean_t["beta"]

library(ggplot2)
library(patchwork)  # install.packages("patchwork") if needed

# Predicted means
mu_female_normal <- mean_normal["alpha"]
mu_female_t <- mean_t["alpha"]

mu_male_normal <- mean_normal["alpha"] + mean_normal["beta"]
mu_male_t <- mean_t["alpha"] + mean_t["beta"]

# --------- Female Plot ---------
p1 <- ggplot(train_data[train_data$sex == "female", ], aes(x = BMI)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "gray90", color = "black") +
  stat_function(aes(color = "Normal"),
                fun = dnorm,
                args = list(mean = mu_female_normal, sd = mean_normal["sigma"]),
                size = 1) +
  stat_function(aes(color = "Student-t"),
                fun = function(x) {
                  dt((x - mu_female_t) / mean_t["sigma"], df = mean_t["nu"]) / mean_t["sigma"]
                },
                size = 1, linetype = "dashed") +
  labs(title = "(a) Female Athletes", x = "BMI", y = "Density", color = "Distribution") +
  theme_minimal(base_size = 13)

# --------- Male Plot ---------
p2 <- ggplot(train_data[train_data$sex == "male", ], aes(x = BMI)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "gray90", color = "black") +
  stat_function(aes(color = "Normal"),
                fun = dnorm,
                args = list(mean = mu_male_normal, sd = mean_normal["sigma"]),
                size = 1) +
  stat_function(aes(color = "Student-t"),
                fun = function(x) {
                  dt((x - mu_male_t) / mean_t["sigma"], df = mean_t["nu"]) / mean_t["sigma"]
                },
                size = 1, linetype = "dashed") +
  labs(title = "(b) Male Athletes", x = "BMI", y = "Density", color = "Distribution") +
  theme_minimal(base_size = 13)

# Combine with shared legend
(p1 | p2) + plot_layout(guides = "collect") & theme(legend.position = "top")





# --------- Prepare Validation Data ---------
x_valid <- valid_data$sex_num
y_valid <- valid_data$BMI
n_valid <- length(y_valid)

# --------- Posterior Samples ---------
# Extract samples as matrices
samp_normal <- as.matrix(samples_normal)
samp_t <- as.matrix(samples_t)
n_samples <- nrow(samp_normal)

# --------- Compute Log Predictive Densities ---------

# Normal model log-likelihoods
log_pred_normal <- matrix(NA, nrow = n_valid, ncol = n_samples)
for (i in 1:n_valid) {
  for (s in 1:n_samples) {
    mu <- samp_normal[s, "alpha"] + samp_normal[s, "beta"] * x_valid[i]
    sd <- samp_normal[s, "sigma"]
    log_pred_normal[i, s] <- dnorm(y_valid[i], mean = mu, sd = sd, log = TRUE)
  }
}

# Student-t model log-likelihoods
log_pred_t <- matrix(NA, nrow = n_valid, ncol = n_samples)
for (i in 1:n_valid) {
  for (s in 1:n_samples) {
    mu <- samp_t[s, "alpha"] + samp_t[s, "beta"] * x_valid[i]
    sigma <- samp_t[s, "sigma"]
    nu <- samp_t[s, "nu"]
    log_pred_t[i, s] <- dt((y_valid[i] - mu) / sigma, df = nu, log = TRUE) - log(sigma)
  }
}

# --------- Compute Average Log Predictive Densities ---------
# Log pointwise predictive density: average over posterior samples
l_normal <- sum(apply(log_pred_normal, 1, mean))
l_t <- sum(apply(log_pred_t, 1, mean))

