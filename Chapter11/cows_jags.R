# =========================================================
# Bayesian Linear Mixed Model: Milk Protein ~ Diet + Time
# Dataset: Milk (nlme package)
# =========================================================

rm(list = ls())

library(nlme)
library(R2jags)
library(coda)
library(dplyr)

# ---------------------------------------------------------
# 1. Data preparation
# ---------------------------------------------------------

data(Milk)
df <- as.data.frame(Milk)
names(df) <- tolower(names(df))

df$diet_lupine   <- ifelse(df$diet == "lupins",        1, 0)
df$diet_barlupin <- ifelse(df$diet == "barley+lupins", 1, 0)
df$cow_id        <- as.integer(factor(df$cow))
df$time_lin      <- df$time
df$time_ch       <- pmax(df$time - 3, 0)
df$time_ch2      <- pmax(df$time - 3, 0)^2

X <- model.matrix(~ diet_lupine + diet_barlupin +
                    time_lin + time_ch + time_ch2,
                  data = df)

data_jags <- list(
  N   = nrow(df),
  M   = length(unique(df$cow_id)),
  K   = ncol(X),
  y   = df$protein,
  X   = X,
  cow = df$cow_id
)

cat("N =", data_jags$N, "\n")
cat("M =", data_jags$M, "\n")
cat("K =", data_jags$K, "\n")

# ---------------------------------------------------------
# 2. JAGS model
# ---------------------------------------------------------

lmm_model <- "
model {

  # Likelihood
  for (n in 1:N) {
    y[n]  ~ dnorm(mu[n], tau_e)
    mu[n] <- inprod(X[n, ], beta[]) + b[cow[n]]
  }

  # Random effects
  for (m in 1:M) {
    b[m] ~ dnorm(0, tau_b)
  }

  # Priors on fixed effects
  for (k in 1:K) {
    beta[k] ~ dnorm(0, 0.001)
  }

  # Priors on precision parameters
  tau_b ~ dgamma(0.001, 0.001)
  tau_e ~ dgamma(0.001, 0.001)

  # Derived quantities
  sigma_b <- 1 / sqrt(tau_b)
  sigma_e <- 1 / sqrt(tau_e)
  icc     <- tau_e / (tau_e + tau_b)
}
"

# ---------------------------------------------------------
# 3. MCMC settings
# ---------------------------------------------------------

params <- c("beta", "sigma_b", "sigma_e", "icc")

inits_fun <- function() {
  list(
    beta  = rnorm(data_jags$K, 0, 0.1),
    b     = rnorm(data_jags$M, 0, 0.1),
    tau_b = rgamma(1, 1, 1),
    tau_e = rgamma(1, 1, 1)
  )
}

# ---------------------------------------------------------
# 4. Model fitting
# ---------------------------------------------------------

set.seed(123)

fit <- jags(
  data               = data_jags,
  inits              = inits_fun,
  parameters.to.save = params,
  model.file         = textConnection(lmm_model),
  n.chains           = 3,
  n.iter             = 10000,
  n.burnin           = 3000,
  n.thin             = 2
)

# ---------------------------------------------------------
# 5. Posterior summary
# ---------------------------------------------------------

print(fit)
round(fit$BUGSoutput$summary, 4)

# ---------------------------------------------------------
# 6. Extract MCMC samples
# ---------------------------------------------------------

sims     <- fit$BUGSoutput$sims.array
n_iter   <- dim(sims)[1]
n_chains <- dim(sims)[2]

# helper: build mcmc.list for one parameter
to_mcmc_list <- function(param) {
  as.mcmc.list(
    lapply(seq_len(n_chains), function(j) mcmc(sims[, j, param]))
  )
}

beta1_mcmc   <- to_mcmc_list("beta[1]")
beta2_mcmc   <- to_mcmc_list("beta[2]")
beta3_mcmc   <- to_mcmc_list("beta[3]")
beta4_mcmc   <- to_mcmc_list("beta[4]")
beta5_mcmc   <- to_mcmc_list("beta[5]")
beta6_mcmc   <- to_mcmc_list("beta[6]")
sigmab_mcmc  <- to_mcmc_list("sigma_b")
sigmae_mcmc  <- to_mcmc_list("sigma_e")
icc_mcmc     <- to_mcmc_list("icc")

# ---------------------------------------------------------
# 7. Trace plots
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

par(mfrow = c(3, 3), mar = c(3, 4, 3, 1))
trace_param(beta1_mcmc,  expression(beta[1]))
trace_param(beta2_mcmc,  expression(beta[2]))
trace_param(beta3_mcmc,  expression(beta[3]))
trace_param(beta4_mcmc,  expression(beta[4]))
trace_param(beta5_mcmc,  expression(beta[5]))
trace_param(beta6_mcmc,  expression(beta[6]))
trace_param(sigmab_mcmc, expression(sigma[b]))
trace_param(sigmae_mcmc, expression(sigma[e]))
trace_param(icc_mcmc,    expression(ICC))
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# ---------------------------------------------------------
# 8. Autocorrelation plots
# ---------------------------------------------------------

acf_list <- list(
  beta1  = autocorr(beta1_mcmc[[1]]),
  beta4  = autocorr(beta4_mcmc[[1]]),
  sigmab = autocorr(sigmab_mcmc[[1]]),
  sigmae = autocorr(sigmae_mcmc[[1]])
)

acf_labels <- list(
  beta1  = expression("ACF -" ~ beta[1]),
  beta4  = expression("ACF -" ~ beta[4]),
  sigmab = expression("ACF -" ~ sigma[b]),
  sigmae = expression("ACF -" ~ sigma[e])
)

lags <- as.numeric(dimnames(acf_list$beta1)[[1]])

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
  abline(h = c(-0.1, 0.1), lty = 2, col = "tomato")
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# ---------------------------------------------------------
# 9. Gelman-Rubin convergence diagnostic
# ---------------------------------------------------------

cat("\n--- Gelman-Rubin diagnostics ---\n")
gr <- lapply(
  list(beta1  = beta1_mcmc,  beta2  = beta2_mcmc,
       beta3  = beta3_mcmc,  beta4  = beta4_mcmc,
       beta5  = beta5_mcmc,  beta6  = beta6_mcmc,
       sigmab = sigmab_mcmc, sigmae = sigmae_mcmc,
       icc    = icc_mcmc),
  gelman.diag
)
lapply(gr, print)

# ---------------------------------------------------------
# 10. Effective sample size
# ---------------------------------------------------------

cat("\n--- Effective sample sizes ---\n")
ess <- sapply(
  list(beta1  = beta1_mcmc,  beta2  = beta2_mcmc,
       beta3  = beta3_mcmc,  beta4  = beta4_mcmc,
       beta5  = beta5_mcmc,  beta6  = beta6_mcmc,
       sigmab = sigmab_mcmc, sigmae = sigmae_mcmc,
       icc    = icc_mcmc),
  effectiveSize
)
print(round(ess))

# ---------------------------------------------------------
# 11. Posterior density plots
# ---------------------------------------------------------

pool <- function(mcmc_list) unlist(lapply(mcmc_list, as.numeric))

post <- list(
  beta1  = pool(beta1_mcmc),
  beta2  = pool(beta2_mcmc),
  beta3  = pool(beta3_mcmc),
  beta4  = pool(beta4_mcmc),
  sigmab = pool(sigmab_mcmc),
  sigmae = pool(sigmae_mcmc)
)

labels <- list(
  beta1  = expression(beta[1] ~ "(intercept)"),
  beta2  = expression(beta[2] ~ "(lupine)"),
  beta3  = expression(beta[3] ~ "(barley+lupine)"),
  beta4  = expression(beta[4] ~ "(time)"),
  sigmab = expression(sigma[b] ~ "(between-cow SD)"),
  sigmae = expression(sigma[e] ~ "(residual SD)")
)

par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (nm in names(post)) {
  d  <- density(post[[nm]])
  ci <- quantile(post[[nm]], c(0.025, 0.975))
  
  plot(d, main = labels[[nm]], xlab = "", lwd = 2,
       col = "steelblue")
  
  x_shade <- c(ci[1], d$x[d$x >= ci[1] & d$x <= ci[2]], ci[2])
  y_shade <- c(0,     d$y[d$x >= ci[1] & d$x <= ci[2]], 0)
  polygon(x_shade, y_shade,
          col = adjustcolor("steelblue", 0.25), border = NA)
  
  abline(v = mean(post[[nm]]), lty = 2, col = "tomato", lwd = 1.5)
  legend("topright",
         legend = c(sprintf("Mean = %.3f", mean(post[[nm]])),
                    sprintf("95%% CI [%.3f, %.3f]", ci[1], ci[2])),
         bty = "n", cex = 0.75)
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# ---------------------------------------------------------
# 12. Clean results table
# ---------------------------------------------------------

beta_names <- c("Intercept", "Lupine", "Barley+Lupine",
                "Time", "Slope change (t>3)", "Quadratic (t>3)")

beta_cols <- paste0("beta[", 1:data_jags$K, "]")
smat      <- fit$BUGSoutput$sims.matrix

results <- data.frame(
  Parameter = beta_names,
  Mean      = round(colMeans(smat[, beta_cols]), 4),
  SD        = round(apply(smat[, beta_cols], 2, sd), 4),
  Lower95   = round(apply(smat[, beta_cols], 2, quantile, 0.025), 4),
  Upper95   = round(apply(smat[, beta_cols], 2, quantile, 0.975), 4),
  row.names = NULL
)
print(results, row.names = FALSE)

cat("\nsigma_b: mean =", round(mean(smat[, "sigma_b"]), 4),
    " 95% CI: [",
    round(quantile(smat[, "sigma_b"], 0.025), 4), ",",
    round(quantile(smat[, "sigma_b"], 0.975), 4), "]\n")

cat("sigma_e: mean =", round(mean(smat[, "sigma_e"]), 4),
    " 95% CI: [",
    round(quantile(smat[, "sigma_e"], 0.025), 4), ",",
    round(quantile(smat[, "sigma_e"], 0.975), 4), "]\n")

cat("ICC:     mean =", round(mean(smat[, "icc"]), 4),
    " 95% CI: [",
    round(quantile(smat[, "icc"], 0.025), 4), ",",
    round(quantile(smat[, "icc"], 0.975), 4), "]\n")

# ---------------------------------------------------------
# 13. DIC
# ---------------------------------------------------------

cat("\n--- DIC ---\n")
cat("DIC =", fit$BUGSoutput$DIC, "\n")
cat("pD  =", fit$BUGSoutput$pD,  "\n")