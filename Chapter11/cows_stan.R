library(rstan)
library(nlme)
library(dplyr)
library(ggplot2)
library(bayesplot)
library(loo)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ── 1. Load and prepare data ───────────────────────────────────────────────
data(Milk)
df <- as.data.frame(Milk)

# rename for convenience
names(df) <- tolower(names(df))   # protein, time, cow, diet

# diet indicators (barley = reference)
df$diet_lupine   <- ifelse(df$diet == "lupine",        1, 0)
df$diet_barlupin <- ifelse(df$diet == "barley+lupine", 1, 0)

# cow index
df$cow_id <- as.integer(factor(df$cow))

# piecewise time terms
df$time_lin <- df$time
df$time_ch  <- pmax(df$time - 3, 0)
df$time_ch2 <- pmax(df$time - 3, 0)^2

cat("N =", nrow(df), "\n")
cat("Cows =", length(unique(df$cow_id)), "\n")
cat("Weeks =", length(unique(df$time)), "\n")
cat("Diet levels:", levels(df$diet), "\n")

 
  
# ── 3. Prepare Stan data ──────────────────────────────────────────────────
X <- model.matrix(~ diet_lupine + diet_barlupin +
                    time_lin + time_ch + time_ch2,
                  data = df)

stan_data <- list(
  N    = nrow(df),
  M    = length(unique(df$cow_id)),
  K    = ncol(X),
  y    = df$protein,
  X    = X,
  cow  = df$cow_id,
  time = df$time
)

# ── 4. Fit Stan model ─────────────────────────────────────────────────────
fit <- stan(
  file    = "cows.stan",
  data    = stan_data,
  chains  = 2,
  iter    = 10000,
  warmup  = 5000,
  seed    = 123,
  control = list(adapt_delta = 0.95)
)

# ── 5. Print posterior summaries ──────────────────────────────────────────
print(fit, pars = c("beta", "sigma_b", "sigma_e", "icc"), digits = 3)

# extract draws
post       <- extract(fit)
beta_mean  <- colMeans(post$beta)
beta_ci    <- apply(post$beta, 2, quantile, probs = c(0.025, 0.975))
beta_names <- c("Intercept", "Lupine", "Barley+Lupine",
                "Time", "Slope change (t>3)", "Quadratic (t>3)")

results <- data.frame(
  Parameter = beta_names,
  Mean      = round(beta_mean, 4),
  SD        = round(apply(post$beta, 2, sd), 4),
  Lower95   = round(beta_ci[1, ], 4),
  Upper95   = round(beta_ci[2, ], 4)
)
print(results, row.names = FALSE)

print(fit, pars = c("beta", "sigma_b", "sigma_e", "icc"))

  
