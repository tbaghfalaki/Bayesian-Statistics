# =========================================================
# Bayesian Multiple Linear Regression: BMI ~ Sex + Blood markers
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
RCC <- ais$RCC
WCC <- ais$WCC
Hc  <- ais$Hc
Hg  <- ais$Hg
Fe  <- ais$Fe

# Remove observations with missing values
complete <- complete.cases(bmi, sex, RCC, WCC, Hc, Hg, Fe)

bmi <- bmi[complete]
sex <- sex[complete]
RCC <- RCC[complete]
WCC <- WCC[complete]
Hc  <- Hc[complete]
Hg  <- Hg[complete]
Fe  <- Fe[complete]

cat("Sample size after removing missing values:", sum(complete), "\n")

data_jags <- list(
  N   = length(bmi),
  bmi = bmi,
  sex = sex,
  RCC = RCC,
  WCC = WCC,
  Hc  = Hc,
  Hg  = Hg,
  Fe  = Fe
)

# ---------------------------------------------------------
# 2. JAGS model
# ---------------------------------------------------------

linear_model <- "
model {

  # Likelihood
  for (i in 1:N) {
    bmi[i] ~ dnorm(mu[i], tau2)
    mu[i]  <- beta0 +
               beta1 * sex[i] +
               beta2 * RCC[i] +
               beta3 * WCC[i] +
               beta4 * Hc[i]  +
               beta5 * Hg[i]  +
               beta6 * Fe[i]
  }

  # Weakly informative priors for regression coefficients
  beta0 ~ dnorm(0, 0.0001)
  beta1 ~ dnorm(0, 0.0001)
  beta2 ~ dnorm(0, 0.0001)
  beta3 ~ dnorm(0, 0.0001)
  beta4 ~ dnorm(0, 0.0001)
  beta5 ~ dnorm(0, 0.0001)
  beta6 ~ dnorm(0, 0.0001)

  # Prior for residual precision
  tau2   ~ dgamma(0.001, 0.001)

  # Derived: residual variance and standard deviation
  sigma2 <- 1 / tau2
}
"

# ---------------------------------------------------------
# 3. MCMC settings
# ---------------------------------------------------------

params <- c(
  "beta0", "beta1", "beta2", "beta3",
  "beta4", "beta5", "beta6",
  "tau2", "sigma2"
)

inits_fun <- function() {
  list(
    beta0 = rnorm(1, mean(bmi), 1),
    beta1 = rnorm(1, 0, 0.5),
    beta2 = rnorm(1, 0, 0.5),
    beta3 = rnorm(1, 0, 0.5),
    beta4 = rnorm(1, 0, 0.5),
    beta5 = rnorm(1, 0, 0.5),
    beta6 = rnorm(1, 0, 0.5),
    tau2  = 1 / var(bmi)
  )
}

N_CHAINS <- 3
N_ITER   <- 15000
N_BURNIN <- 3000
N_THIN   <- 3

# ---------------------------------------------------------
# 4. Model fitting
# ---------------------------------------------------------

set.seed(42)

fit <- jags(
  data               = data_jags,
  inits              = inits_fun,
  parameters.to.save = params,
  model.file         = textConnection(linear_model),
  n.chains           = N_CHAINS,
  n.iter             = N_ITER,
  n.burnin           = N_BURNIN,
  n.thin             = N_THIN
)

# ---------------------------------------------------------
# 5. Posterior summary
# ---------------------------------------------------------

print(fit)
round(fit$BUGSoutput$summary, 4)

# ---------------------------------------------------------
# 6. Extract MCMC samples
# ---------------------------------------------------------

sims     <- fit$BUGSoutput$sims.array   # [iterations x chains x parameters]
n_iter   <- dim(sims)[1]
n_chains <- dim(sims)[2]

# Helper: build mcmc.list for one parameter
to_mcmc_list <- function(param) {
  as.mcmc.list(
    lapply(seq_len(n_chains), function(j) mcmc(sims[, j, param]))
  )
}

beta0_mcmc  <- to_mcmc_list("beta0")
beta1_mcmc  <- to_mcmc_list("beta1")
beta2_mcmc  <- to_mcmc_list("beta2")
beta3_mcmc  <- to_mcmc_list("beta3")
beta4_mcmc  <- to_mcmc_list("beta4")
beta5_mcmc  <- to_mcmc_list("beta5")
beta6_mcmc  <- to_mcmc_list("beta6")
tau2_mcmc   <- to_mcmc_list("tau2")
sigma2_mcmc <- to_mcmc_list("sigma2")

all_mcmc <- list(
  beta0  = beta0_mcmc,
  beta1  = beta1_mcmc,
  beta2  = beta2_mcmc,
  beta3  = beta3_mcmc,
  beta4  = beta4_mcmc,
  beta5  = beta5_mcmc,
  beta6  = beta6_mcmc,
  tau2   = tau2_mcmc,
  sigma2 = sigma2_mcmc
)

param_labels <- list(
  beta0  = expression(beta[0] ~ "(intercept)"),
  beta1  = expression(beta[1] ~ "(sex)"),
  beta2  = expression(beta[2] ~ "(RCC)"),
  beta3  = expression(beta[3] ~ "(WCC)"),
  beta4  = expression(beta[4] ~ "(Hc)"),
  beta5  = expression(beta[5] ~ "(Hg)"),
  beta6  = expression(beta[6] ~ "(Fe)"),
  tau2   = expression(tau^2),
  sigma2 = expression(sigma^2)
)

param_names <- names(all_mcmc)   # 9 parameters

# ---------------------------------------------------------
# 7. Trace plots  (all chains overlaid, all 9 in one 5x2 grid)
# ---------------------------------------------------------

chain_cols <- c("steelblue", "tomato", "forestgreen")

trace_param <- function(mcmc_list, param_label) {
  plot(
    NA,
    xlim = c(1, n_iter),
    ylim = range(unlist(mcmc_list)),
    type = "n",
    main = param_label,
    xlab = "Iteration",
    ylab = ""
  )
  for (j in seq_len(n_chains))
    lines(as.numeric(mcmc_list[[j]]), col = chain_cols[j], lwd = 0.6)
}

# All 9 parameters in one 3x3 panel
par(mfrow = c(3, 3), mar = c(3, 3, 2, 1), oma = c(0, 0, 3, 0))
for (nm in param_names) {
  trace_param(all_mcmc[[nm]], param_labels[[nm]])
}
# Shared title and legend in outer margin
mtext("Trace plots (all chains)", outer = TRUE, cex = 1.2, font = 2)
legend("topright", legend = paste("Chain", seq_len(n_chains)),
       col = chain_cols, lty = 1, bty = "n", cex = 0.9, xpd = NA)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0))

# ---------------------------------------------------------
# 8. Autocorrelation plots - chain 1 only (2x2 per page)
# ---------------------------------------------------------

# autocorr() returns values without overriding par()
acf_vals <- lapply(all_mcmc, function(m) autocorr(m[[1]]))
lags      <- as.numeric(dimnames(acf_vals[[1]])[[1]])

# All 9 parameters in one 3x3 panel
par(mfrow = c(3, 3), mar = c(4, 4, 2, 1), oma = c(0, 0, 3, 0))
for (nm in param_names) {
  av <- acf_vals[[nm]][, 1, 1]
  barplot(
    av,
    names.arg = lags,
    main      = param_labels[[nm]],
    xlab      = "Lag",
    ylab      = "ACF",
    col       = "steelblue",
    border    = NA,
    ylim      = c(-1, 1)
  )
  abline(h = 0,             lty = 1)
  abline(h = c(-0.1, 0.1), lty = 2, col = "tomato")
}
mtext("Autocorrelation plots - Chain 1", outer = TRUE, cex = 1.2, font = 2)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0))

# ---------------------------------------------------------
# 9. Gelman-Rubin convergence diagnostic  (R-hat < 1.1 is good)
# ---------------------------------------------------------

cat("\n--- Gelman-Rubin diagnostics ---\n")
gr <- lapply(all_mcmc, gelman.diag)
lapply(gr, print)

# ---------------------------------------------------------
# 10. Effective sample size
# ---------------------------------------------------------

cat("\n--- Effective sample sizes ---\n")
ess <- sapply(all_mcmc, effectiveSize)
print(round(ess))

# ---------------------------------------------------------
# 11. Posterior density plots (2x2 per page, with 95% CI)
# ---------------------------------------------------------

pool <- function(mcmc_list) unlist(lapply(mcmc_list, as.numeric))
post <- lapply(all_mcmc, pool)

# All 9 parameters in one 3x3 panel
par(mfrow = c(3, 3), mar = c(4, 4, 2, 1), oma = c(0, 0, 3, 0))
for (nm in param_names) {
  d  <- density(post[[nm]])
  ci <- quantile(post[[nm]], c(0.025, 0.975))
  
  plot(d, main = param_labels[[nm]], xlab = "", lwd = 2, col = "steelblue")
  
  x_shade <- c(ci[1], d$x[d$x >= ci[1] & d$x <= ci[2]], ci[2])
  y_shade <- c(0,     d$y[d$x >= ci[1] & d$x <= ci[2]], 0)
  polygon(x_shade, y_shade, col = adjustcolor("steelblue", 0.25), border = NA)
  
  abline(v = mean(post[[nm]]), lty = 2, col = "tomato", lwd = 1.5)
  legend("topright",
         legend = c(sprintf("Mean = %.3f", mean(post[[nm]])),
                    sprintf("95%% CI [%.3f, %.3f]", ci[1], ci[2])),
         bty = "n", cex = 0.7)
}
mtext("Posterior density plots (95% CI shaded)", outer = TRUE, cex = 1.2, font = 2)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0))