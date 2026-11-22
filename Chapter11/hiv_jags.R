# =========================================================
# Bayesian Linear Mixed Model: sqrt(CD4) ~ Drug + Time
# Dataset: aids (JM package) — Goldman et al. (1996)
# JAGS implementation
# =========================================================

rm(list = ls())

library(JM)
library(R2jags)
library(coda)
library(ggplot2)
library(dplyr)

# ---------------------------------------------------------
# 1. Load and prepare data
# ---------------------------------------------------------

data(aids)
df <- aids

cat("Patients:", length(unique(df$patient)), "\n")
cat("Obs:",      nrow(df),                   "\n")
cat("Times:",    unique(df$obstime),          "\n")
table(df$drug)
table(df$gender)
table(df$prevOI)
table(df$AZT)

df$sqrtCD4    <- sqrt(df$CD4)
df$patient_id <- as.integer(factor(df$patient))
df$drug_i     <- as.integer(df$drug    == "ddI")
df$gender_i   <- as.integer(df$gender  == "male")
df$prevOI_i   <- as.integer(df$prevOI  == "AIDS")
df$AZT_i      <- as.integer(df$AZT     == "failure")

# ---------------------------------------------------------
# 2. Spaghetti plot
# ---------------------------------------------------------

drug_colours <- c("ddC" = "#2166AC",
                  "ddI" = "#D6604D")

mean_traj <- df %>%
  group_by(drug, obstime) %>%
  summarise(mean_cd4 = mean(sqrtCD4), .groups = "drop")

p_spaghetti <- ggplot(df, aes(x      = obstime,
                              y      = sqrtCD4,
                              group  = patient,
                              colour = drug)) +
  geom_line(alpha = 0.25, linewidth = 0.35) +
  geom_line(data        = mean_traj,
            aes(x       = obstime,
                y       = mean_cd4,
                group   = drug,
                colour  = drug),
            linewidth   = 1.6,
            alpha       = 1.0,
            inherit.aes = FALSE) +
  geom_vline(xintercept = c(2, 6, 12),
             linetype   = "dashed",
             colour     = "grey50",
             linewidth  = 0.4) +
  facet_wrap(~ drug, ncol = 2) +
  scale_colour_manual(values = drug_colours) +
  scale_x_continuous(breaks = c(0, 2, 6, 12, 18)) +
  labs(x      = "Observation time (months)",
       y      = expression(sqrt(CD4)),
       colour = "Drug") +
  theme_bw(base_size = 11) +
  theme(legend.position  = "none",
        strip.background = element_rect(fill = "grey92"),
        strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank())

print(p_spaghetti)

# ---------------------------------------------------------
# 3. JAGS data list
# ---------------------------------------------------------

data_jags <- list(
  N       = nrow(df),
  M       = length(unique(df$patient_id)),
  y       = df$sqrtCD4,
  time    = df$obstime,
  drug    = df$drug_i,
  gender  = df$gender_i,
  prevOI  = df$prevOI_i,
  AZT     = df$AZT_i,
  patient = df$patient_id,
  R       = diag(2)              # scale matrix for Wishart prior
)

cat("N =", data_jags$N, "\n")
cat("M =", data_jags$M, "\n")

# ---------------------------------------------------------
# 4. JAGS model
# ---------------------------------------------------------

lmm_aids <- "
model {

  # ── likelihood ───────────────────────────────────────
  for (n in 1:N) {
    y[n]  ~ dnorm(mu[n], tau_e)
    mu[n] <- beta[1]
           + beta[2] * time[n]
           + beta[3] * time[n] * drug[n]
           + beta[4] * gender[n]
           + beta[5] * prevOI[n]
           + beta[6] * AZT[n]
           + b[patient[n], 1]
           + b[patient[n], 2] * time[n]
  }

  # ── random effects ───────────────────────────────────
  for (m in 1:M) {
    b[m, 1:2] ~ dmnorm(mu0[], Omega[,])
  }

  # mean of random effects (zero)
  mu0[1] <- 0
  mu0[2] <- 0

  # ── priors on fixed effects ──────────────────────────
  for (k in 1:6) {
    beta[k] ~ dnorm(0, 0.01)
  }

  # ── prior on precision matrix (Wishart) ──────────────
  Omega[1:2, 1:2] ~ dwish(R[,], 3)

  # ── prior on residual precision ──────────────────────
  tau_e ~ dgamma(0.001, 0.001)

  # ── derived quantities ───────────────────────────────
  sigma_e    <- 1 / sqrt(tau_e)
  D[1:2,1:2] <- inverse(Omega[,])
  sigma_b0   <- sqrt(D[1, 1])
  sigma_b1   <- sqrt(D[2, 2])
  rho        <- D[1, 2] / (sigma_b0 * sigma_b1)
  icc        <- D[1, 1] / (D[1, 1] + 1 / tau_e)
}
"

writeLines(lmm_aids, "aids_lmm.jags")

# ---------------------------------------------------------
# 5. Initial values
# ---------------------------------------------------------

inits_fun <- function() {
  list(
    beta  = rnorm(6, 0, 0.1),
    tau_e = 1,
    Omega = diag(2),
    b     = matrix(rnorm(data_jags$M * 2, 0, 0.1),
                   nrow = data_jags$M, ncol = 2)
  )
}

# ---------------------------------------------------------
# 6. Fit model
# ---------------------------------------------------------

params <- c("beta", "sigma_b0", "sigma_b1",
            "sigma_e", "rho", "icc", "D")

set.seed(123)

fit <- jags(
  data               = data_jags,
  inits              = inits_fun,
  parameters.to.save = params,
  model.file         = textConnection(lmm_aids),
  n.chains           = 3,
  n.iter             = 10000,
  n.burnin           = 3000,
  n.thin             = 2
)

# ---------------------------------------------------------
# 7. Posterior summaries
# ---------------------------------------------------------

print(fit)
round(fit$BUGSoutput$summary[
  c("beta[1]","beta[2]","beta[3]",
    "beta[4]","beta[5]","beta[6]",
    "sigma_b0","sigma_b1",
    "sigma_e","rho","icc"), ], 3)

# ---------------------------------------------------------
# 8. Extract MCMC samples
# ---------------------------------------------------------

sims     <- fit$BUGSoutput$sims.array
n_iter   <- dim(sims)[1]
n_chains <- dim(sims)[2]

to_mcmc_list <- function(param) {
  as.mcmc.list(
    lapply(seq_len(n_chains),
           function(j) mcmc(sims[, j, param]))
  )
}

beta1_mcmc  <- to_mcmc_list("beta[1]")
beta2_mcmc  <- to_mcmc_list("beta[2]")
beta5_mcmc  <- to_mcmc_list("beta[5]")
sigmab_mcmc <- to_mcmc_list("sigma_b0")
sigmae_mcmc <- to_mcmc_list("sigma_e")
rho_mcmc    <- to_mcmc_list("rho")

# ---------------------------------------------------------
# 9. Trace plots
# ---------------------------------------------------------

chain_cols <- c("steelblue", "tomato", "forestgreen")

trace_param <- function(mcmc_list, param_label) {
  plot(NA,
       xlim = c(1, n_iter),
       ylim = range(unlist(mcmc_list)),
       type = "n",
       main = bquote("Trace -" ~ .(param_label)),
       xlab = "Iteration",
       ylab = param_label)
  for (j in seq_len(n_chains))
    lines(as.numeric(mcmc_list[[j]]),
          col = chain_cols[j], lwd = 0.6)
  legend("topright",
         legend = paste("Chain", seq_len(n_chains)),
         col    = chain_cols,
         lty    = 1, bty = "n", cex = 0.7)
}

par(mfrow = c(2, 3), mar = c(3, 4, 3, 1))
trace_param(beta1_mcmc,  expression(beta[1]))
trace_param(beta2_mcmc,  expression(beta[2]))
trace_param(beta5_mcmc,  expression(beta[5]))
trace_param(sigmab_mcmc, expression(sigma[b[0]]))
trace_param(sigmae_mcmc, expression(sigma[e]))
trace_param(rho_mcmc,    expression(rho))
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# ---------------------------------------------------------
# 10. Autocorrelation plots
# ---------------------------------------------------------

acf_list <- list(
  beta1  = autocorr(beta1_mcmc[[1]]),
  beta2  = autocorr(beta2_mcmc[[1]]),
  sigmab = autocorr(sigmab_mcmc[[1]]),
  rho    = autocorr(rho_mcmc[[1]])
)

acf_labels <- list(
  beta1  = expression("ACF -" ~ beta[1]),
  beta2  = expression("ACF -" ~ beta[2]),
  sigmab = expression("ACF -" ~ sigma[b[0]]),
  rho    = expression("ACF -" ~ rho)
)

lags <- as.numeric(dimnames(acf_list$beta1)[[1]])

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (nm in names(acf_list)) {
  acf_vals <- acf_list[[nm]][, 1, 1]
  barplot(acf_vals,
          names.arg = lags,
          main      = acf_labels[[nm]],
          xlab      = "Lag",
          ylab      = "ACF",
          col       = "steelblue",
          border    = NA,
          ylim      = c(-1, 1))
  abline(h = 0, lty = 1)
  abline(h = c(-0.1, 0.1), lty = 2, col = "tomato")
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# ---------------------------------------------------------
# 11. Gelman-Rubin diagnostics
# ---------------------------------------------------------

cat("\n--- Gelman-Rubin diagnostics ---\n")
gr <- lapply(
  list(beta1  = beta1_mcmc,  beta2  = beta2_mcmc,
       beta5  = beta5_mcmc,  sigmab = sigmab_mcmc,
       sigmae = sigmae_mcmc, rho    = rho_mcmc),
  gelman.diag
)
lapply(gr, print)

# ---------------------------------------------------------
# 12. Effective sample sizes
# ---------------------------------------------------------

cat("\n--- Effective sample sizes ---\n")
ess <- sapply(
  list(beta1  = beta1_mcmc,  beta2  = beta2_mcmc,
       beta5  = beta5_mcmc,  sigmab = sigmab_mcmc,
       sigmae = sigmae_mcmc, rho    = rho_mcmc),
  effectiveSize
)
print(round(ess))

# ---------------------------------------------------------
# 13. Posterior density plots
# ---------------------------------------------------------

pool <- function(mcmc_list) unlist(lapply(mcmc_list, as.numeric))

post_pool <- list(
  beta1  = pool(beta1_mcmc),
  beta2  = pool(beta2_mcmc),
  beta5  = pool(beta5_mcmc),
  sigmab = pool(sigmab_mcmc),
  sigmae = pool(sigmae_mcmc),
  rho    = pool(rho_mcmc)
)

dens_labels <- list(
  beta1  = expression(beta[1] ~ "(intercept)"),
  beta2  = expression(beta[2] ~ "(time)"),
  beta5  = expression(beta[5] ~ "(prior OI)"),
  sigmab = expression(sigma[b[0]] ~ "(intercept SD)"),
  sigmae = expression(sigma[e] ~ "(residual SD)"),
  rho    = expression(rho ~ "(intercept-slope corr.)")
)

par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (nm in names(post_pool)) {
  d  <- density(post_pool[[nm]])
  ci <- quantile(post_pool[[nm]], c(0.025, 0.975))
  plot(d, main = dens_labels[[nm]],
       xlab = "", lwd = 2, col = "steelblue")
  x_shade <- c(ci[1],
               d$x[d$x >= ci[1] & d$x <= ci[2]],
               ci[2])
  y_shade <- c(0,
               d$y[d$x >= ci[1] & d$x <= ci[2]],
               0)
  polygon(x_shade, y_shade,
          col    = adjustcolor("steelblue", 0.25),
          border = NA)
  abline(v   = mean(post_pool[[nm]]),
         lty = 2, col = "tomato", lwd = 1.5)
  legend("topright",
         legend = c(sprintf("Mean = %.3f",
                            mean(post_pool[[nm]])),
                    sprintf("95%% CI [%.3f, %.3f]",
                            ci[1], ci[2])),
         bty = "n", cex = 0.75)
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# ---------------------------------------------------------
# 14. Clean results table
# ---------------------------------------------------------

smat <- fit$BUGSoutput$sims.matrix

beta_names <- c("Intercept", "Time",
                "Drug x Time", "Gender",
                "Prior OI", "AZT")

beta_cols <- paste0("beta[", 1:6, "]")

results <- data.frame(
  Parameter = c(beta_names,
                "sigma_b0", "sigma_b1",
                "sigma_e",  "rho", "icc"),
  Mean    = round(colMeans(
    smat[, c(beta_cols,
             "sigma_b0","sigma_b1",
             "sigma_e","rho","icc")]), 3),
  SD      = round(apply(
    smat[, c(beta_cols,
             "sigma_b0","sigma_b1",
             "sigma_e","rho","icc")], 2, sd), 3),
  Lower95 = round(apply(
    smat[, c(beta_cols,
             "sigma_b0","sigma_b1",
             "sigma_e","rho","icc")],
    2, quantile, 0.025), 3),
  Upper95 = round(apply(
    smat[, c(beta_cols,
             "sigma_b0","sigma_b1",
             "sigma_e","rho","icc")],
    2, quantile, 0.975), 3),
  row.names = NULL
)
print(results, row.names = FALSE)

# ---------------------------------------------------------
# 15. DIC
# ---------------------------------------------------------

cat("\n--- DIC ---\n")
cat("DIC =", fit$BUGSoutput$DIC, "\n")
cat("pD  =", fit$BUGSoutput$pD,  "\n")