# =========================================================
# Bayesian Linear Regression: BMI ~ Sex
# Dataset: AIS (Australian Institute of Sport)
# =========================================================

rm(list = ls())

library(sn)
library(R2jags)
library(coda)

# ---------------------------------------------------------
# 1. Data preparation
# ---------------------------------------------------------

data(ais)

bmi <- ais$BMI
sex <- ifelse(ais$sex == "female", 0, 1)   # 0 = female, 1 = male

data_jags <- list(
  N   = length(bmi),
  bmi = bmi,
  sex = sex
)

# ---------------------------------------------------------
# 2. JAGS model
# ---------------------------------------------------------

linear_model <- "
model {

  # Likelihood
  for (i in 1:N) {
    bmi[i] ~ dnorm(mu[i], tau2)
    mu[i]  <- beta0 + beta1 * sex[i]
  }

  # Priors (weakly informative)
  beta0 ~ dnorm(0, 0.0001)
  beta1 ~ dnorm(0, 0.0001)

  # Precision & derived standard deviation
  tau2  ~ dgamma(0.001, 0.001)
  sigma2 <- 1 / tau2
}
"

# ---------------------------------------------------------
# 3. MCMC settings
# ---------------------------------------------------------

params <- c("beta0", "beta1", "tau2", "sigma2")

inits_fun <- function() {
  list(
    beta0 = rnorm(1, mean(bmi), 1),   # randomised for better chain diversity
    beta1 = rnorm(1, 0, 0.5),
    tau2  = 1 / var(bmi)
  )
}


# ---------------------------------------------------------
# 4. Model fitting
# ---------------------------------------------------------

set.seed(42)

fit <- jags(
  data                = data_jags,
  inits               = inits_fun,
  parameters.to.save  = params,
  model.file          = textConnection(linear_model),
  n.chains            = 3,
  n.iter              = 10000,
  n.burnin            = 3000,
  n.thin              = 2
)

# ---------------------------------------------------------
# 5. Posterior summary
# ---------------------------------------------------------

print(fit)
round(fit$BUGSoutput$summary, 4)

# ---------------------------------------------------------
# 6. Extract MCMC samples
# ---------------------------------------------------------

sims <- fit$BUGSoutput$sims.array   # [iterations x chains x parameters]
n_iter   <- dim(sims)[1]
n_chains <- dim(sims)[2]

# Helper: build mcmc.list for one parameter
to_mcmc_list <- function(param) {
  as.mcmc.list(
    lapply(seq_len(n_chains), function(j) mcmc(sims[, j, param]))
  )
}

beta0_mcmc <- to_mcmc_list("beta0")
beta1_mcmc <- to_mcmc_list("beta1")
tau2_mcmc  <- to_mcmc_list("tau2")
sigma2_mcmc <- to_mcmc_list("sigma2")

# ---------------------------------------------------------
# 7. Trace plots  (all chains overlaid, fixed bug)
# ---------------------------------------------------------

chain_cols <- c("steelblue", "tomato", "forestgreen")

trace_param <- function(mcmc_list, param_label) {
  plot(
    NA,
    xlim = c(1, n_iter),
    ylim = range(unlist(mcmc_list)),
    type = "n",
    main = bquote("Trace plot -" ~ .(param_label)),
    xlab = "Iteration",
    ylab = param_label
  )
  for (j in seq_len(n_chains))
    lines(as.numeric(mcmc_list[[j]]), col = chain_cols[j], lwd = 0.6)
  legend("topright", legend = paste("Chain", seq_len(n_chains)),
         col = chain_cols, lty = 1, bty = "n", cex = 0.75)
}

par(mfrow = c(4, 1), mar = c(3, 4, 3, 1))
trace_param(beta0_mcmc, expression(beta[0]))
trace_param(beta1_mcmc, expression(beta[1]))
trace_param(tau2_mcmc,  expression(tau^2))
trace_param(sigma2_mcmc, expression(sigma^2))
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# ---------------------------------------------------------
# 8. Autocorrelation plots
# ---------------------------------------------------------
# ---------------------------------------------------------
# 8. Autocorrelation plots - chain 1 only (2x2 layout)
# ---------------------------------------------------------

acf_list <- list(
  beta0 = autocorr(beta0_mcmc[[1]]),
  beta1 = autocorr(beta1_mcmc[[1]]),
  tau2  = autocorr(tau2_mcmc[[1]]),
  sigma2 = autocorr(sigma2_mcmc[[1]])
)

acf_labels <- list(
  beta0 = expression("ACF -" ~ beta[0]),
  beta1 = expression("ACF -" ~ beta[1]),
  tau2  = expression("ACF -" ~ tau^2),
  sigma2 = expression("ACF -" ~ sigma^2)
)

lags <- as.numeric(dimnames(acf_list$beta0)[[1]])

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (nm in names(acf_list)) {
  acf_vals <- acf_list[[nm]][, 1, 1]
  barplot(
    acf_vals,
    names.arg = lags,
    main      = acf_labels[[nm]],
    xlab      = "Lag",
    ylab      = "ACF",
    col       = "steelblue",
    border    = NA,
    ylim      = c(-1, 1)
  )
  abline(h = 0, lty = 1)
  abline(h = c(-0.1, 0.1), lty = 2, col = "tomato")  # rough ±0.1 reference lines
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# ---------------------------------------------------------
# 9. Gelman-Rubin convergence diagnostic  (R-hat < 1.1 is good)
# ---------------------------------------------------------

cat("\n--- Gelman-Rubin diagnostics ---\n")
gr <- lapply(
  list(beta0 = beta0_mcmc, beta1 = beta1_mcmc,
       tau2  = tau2_mcmc,  sigma2 = sigma2_mcmc),
  gelman.diag
)
lapply(gr, print)

# ---------------------------------------------------------
# 10. Effective sample size
# ---------------------------------------------------------

cat("\n--- Effective sample sizes ---\n")
ess <- sapply(
  list(beta0 = beta0_mcmc, beta1 = beta1_mcmc,
       tau2  = tau2_mcmc,  sigma = sigma2_mcmc),
  effectiveSize
)
print(round(ess))

# ---------------------------------------------------------
# 11. Posterior density plots
# ---------------------------------------------------------

# Pool all chains into one vector per parameter
pool <- function(mcmc_list) unlist(lapply(mcmc_list, as.numeric))

post <- list(
  beta0 = pool(beta0_mcmc),
  beta1 = pool(beta1_mcmc),
  tau2  = pool(tau2_mcmc),
  sigma = pool(sigma2_mcmc)
)

labels <- list(
  beta0 = expression(beta[0] ~ "(intercept)"),
  beta1 = expression(beta[1] ~ "(sex effect)"),
  tau2  = expression(tau^2 ~ "(precision)"),
  sigma = expression(sigma^2 ~ "(variance)")
)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (nm in names(post)) {
  d <- density(post[[nm]])
  plot(d, main = labels[[nm]], xlab = "", lwd = 2, col = "steelblue")
  
  # 95% credible interval shading
  ci <- quantile(post[[nm]], c(0.025, 0.975))
  x_shade <- c(ci[1], d$x[d$x >= ci[1] & d$x <= ci[2]], ci[2])
  y_shade <- c(0,     d$y[d$x >= ci[1] & d$x <= ci[2]], 0)
  polygon(x_shade, y_shade, col = adjustcolor("steelblue", 0.25), border = NA)
  
  abline(v = mean(post[[nm]]), lty = 2, col = "tomato", lwd = 1.5)
  legend("topright",
         legend = c(sprintf("Mean = %.3f", mean(post[[nm]])),
                    sprintf("95%% CI [%.3f, %.3f]", ci[1], ci[2])),
         bty = "n", cex = 0.75)
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
