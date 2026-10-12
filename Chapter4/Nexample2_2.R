# ============================================================
# 1. Load packages and AIS data
# ============================================================

library(sn)
library(R2jags)
library(ggplot2)
library(patchwork)

data(ais, package = "sn")

# Encode sex: 0 = female, 1 = male
ais$sex_num <- ifelse(ais$sex == "male", 1, 0)

# Split into training (70%) and validation (30%) sets
set.seed(123)
n <- nrow(ais)

train_indices <- sample(
  seq_len(n),
  size = round(0.7 * n)
)

train_data <- ais[train_indices, ]
valid_data <- ais[-train_indices, ]

# Data for JAGS
jags_data <- list(
  y = train_data$BMI,
  x = train_data$sex_num,
  N = nrow(train_data)
)

# ============================================================
# 2. Normal regression model
# ============================================================

model_normal <- "
model {
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + beta * x[i]
  }

  alpha ~ dnorm(0, 0.001)
  beta  ~ dnorm(0, 0.001)

  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 100)
}
"

# R2jags requires a model file
file_normal <- tempfile(fileext = ".bug")
writeLines(model_normal, file_normal)

set.seed(123)

fit_normal <- R2jags::jags(
  data = jags_data,
  parameters.to.save = c("alpha", "beta", "sigma"),
  model.file = file_normal,
  n.chains = 3,
  n.iter = 6000,
  n.burnin = 1000,
  n.thin = 1,
  DIC = FALSE
)

# Posterior summaries
print(fit_normal$BUGSoutput$summary)

# Extract all retained posterior draws
samp_normal <- fit_normal$BUGSoutput$sims.matrix

# ============================================================
# 3. Student-t regression model
# ============================================================

model_t <- "
model {
  for (i in 1:N) {
    y[i] ~ dt(mu[i], tau, nu)
    mu[i] <- alpha + beta * x[i]
  }

  alpha ~ dnorm(0, 0.001)
  beta  ~ dnorm(0, 0.001)

  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 100)
  nu ~ dunif(1, 30)
}
"

file_t <- tempfile(fileext = ".bug")
writeLines(model_t, file_t)

set.seed(456)

fit_t <- R2jags::jags(
  data = jags_data,
  parameters.to.save = c("alpha", "beta", "sigma", "nu"),
  model.file = file_t,
  n.chains = 3,
  n.iter = 6000,
  n.burnin = 1000,
  n.thin = 1,
  DIC = FALSE
)

# Posterior summaries
print(fit_t$BUGSoutput$summary)

# Extract all retained posterior draws
samp_t <- fit_t$BUGSoutput$sims.matrix

# Remove temporary model files
unlink(c(file_normal, file_t))

# ============================================================
# 4. Fitted densities by sex
#    Evaluated at posterior mean parameter estimates
# ============================================================

mean_normal <- colMeans(samp_normal)
mean_t <- colMeans(samp_t)

# Female: x = 0
mu_female_normal <- mean_normal["alpha"]
mu_female_t <- mean_t["alpha"]

# Male: x = 1
mu_male_normal <- mean_normal["alpha"] + mean_normal["beta"]
mu_male_t <- mean_t["alpha"] + mean_t["beta"]

# Female plot
p1 <- ggplot(
  train_data[train_data$sex == "female", ],
  aes(x = BMI)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = "gray90",
    color = "black"
  ) +
  stat_function(
    aes(color = "Normal"),
    fun = dnorm,
    args = list(
      mean = mu_female_normal,
      sd = mean_normal["sigma"]
    ),
    linewidth = 1
  ) +
  stat_function(
    aes(color = "Student-t"),
    fun = function(x) {
      dt(
        (x - mu_female_t) / mean_t["sigma"],
        df = mean_t["nu"]
      ) / mean_t["sigma"]
    },
    linewidth = 1,
    linetype = "dashed"
  ) +
  labs(
    title = "(a) Female Athletes",
    x = "BMI",
    y = "Density",
    color = "Distribution"
  ) +
  theme_minimal(base_size = 13)

# Male plot
p2 <- ggplot(
  train_data[train_data$sex == "male", ],
  aes(x = BMI)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = "gray90",
    color = "black"
  ) +
  stat_function(
    aes(color = "Normal"),
    fun = dnorm,
    args = list(
      mean = mu_male_normal,
      sd = mean_normal["sigma"]
    ),
    linewidth = 1
  ) +
  stat_function(
    aes(color = "Student-t"),
    fun = function(x) {
      dt(
        (x - mu_male_t) / mean_t["sigma"],
        df = mean_t["nu"]
      ) / mean_t["sigma"]
    },
    linewidth = 1,
    linetype = "dashed"
  ) +
  labs(
    title = "(b) Male Athletes",
    x = "BMI",
    y = "Density",
    color = "Distribution"
  ) +
  theme_minimal(base_size = 13)

# Combine plots with a shared legend
(p1 | p2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

# ============================================================
# 5. Log posterior predictive density on validation data
# ============================================================

x_valid <- valid_data$sex_num
y_valid <- valid_data$BMI
n_valid <- length(y_valid)

# Numerically stable log-mean-exp function
log_mean_exp <- function(z) {
  m <- max(z)
  m + log(mean(exp(z - m)))
}

# Pointwise log posterior predictive densities
lpd_normal <- numeric(n_valid)
lpd_t <- numeric(n_valid)

for (i in seq_len(n_valid)) {
  
  # ----------------------------------------------------------
  # Normal model
  # ----------------------------------------------------------
  
  mu_normal <- samp_normal[, "alpha"] +
    samp_normal[, "beta"] * x_valid[i]
  
  log_dens_normal <- dnorm(
    y_valid[i],
    mean = mu_normal,
    sd = samp_normal[, "sigma"],
    log = TRUE
  )
  
  lpd_normal[i] <- log_mean_exp(log_dens_normal)
  
  # ----------------------------------------------------------
  # Student-t model
  # ----------------------------------------------------------
  
  mu_t <- samp_t[, "alpha"] +
    samp_t[, "beta"] * x_valid[i]
  
  log_dens_t <- dt(
    (y_valid[i] - mu_t) / samp_t[, "sigma"],
    df = samp_t[, "nu"],
    log = TRUE
  ) - log(samp_t[, "sigma"])
  
  lpd_t[i] <- log_mean_exp(log_dens_t)
}

# ============================================================
# 6. Compare predictive performance
# ============================================================

results <- data.frame(
  Model = c("Normal", "Student-t"),
  Total_LPPD = c(
    sum(lpd_normal),
    sum(lpd_t)
  ),
  Mean_LPPD = c(
    mean(lpd_normal),
    mean(lpd_t)
  )
)

print(results, row.names = FALSE)


