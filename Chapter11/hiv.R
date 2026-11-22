rm(list = ls())

library(JM)
library(rstan)
library(loo)
library(ggplot2)
library(dplyr)
library(bayesplot)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

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

# sqrt transformation
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

p_spaghetti <- ggplot(df, aes(x     = obstime,
                              y     = sqrtCD4,
                              group = patient,
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
# 3. Stan data list
# ---------------------------------------------------------

stan_data <- list(
  N       = nrow(df),
  M       = length(unique(df$patient_id)),
  y       = df$sqrtCD4,
  time    = df$obstime,
  drug    = df$drug_i,
  gender  = df$gender_i,
  prevOI  = df$prevOI_i,
  AZT     = df$AZT_i,
  patient = df$patient_id
)

cat("N =", stan_data$N, "\n")
cat("M =", stan_data$M, "\n")

# ---------------------------------------------------------
# 4. Stan model (save as aids_lmm.stan)
# ---------------------------------------------------------

stan_model_code <- "
data {
  int<lower=1> N;
  int<lower=1> M;
  vector[N] y;
  vector[N] time;
  array[N] int<lower=0, upper=1> drug;
  array[N] int<lower=0, upper=1> gender;
  array[N] int<lower=0, upper=1> prevOI;
  array[N] int<lower=0, upper=1> AZT;
  array[N] int<lower=1, upper=M> patient;
}

parameters {
  real beta1;
  real beta2;
  real beta3;
  real beta4;
  real beta5;
  real beta6;
  matrix[M, 2] b;
  cholesky_factor_corr[2] L_corr;
  vector<lower=0>[2] sigma_b;
  real<lower=0> sigma_e;
}

transformed parameters {
  matrix[2, 2] D;
  D = diag_pre_multiply(sigma_b, L_corr) *
      diag_pre_multiply(sigma_b, L_corr)';
  vector[N] mu;
  for (n in 1:N) {
    mu[n] = beta1
          + beta2  * time[n]
          + beta3  * (time[n] * drug[n])
          + beta4  * gender[n]
          + beta5  * prevOI[n]
          + beta6  * AZT[n]
          + b[patient[n], 1]
          + b[patient[n], 2] * time[n];
  }
}

model {
  beta1 ~ normal(0, 10);
  beta2 ~ normal(0, 10);
  beta3 ~ normal(0, 10);
  beta4 ~ normal(0, 10);
  beta5 ~ normal(0, 10);
  beta6 ~ normal(0, 10);
  sigma_b ~ cauchy(0, 2.5);
  sigma_e ~ cauchy(0, 2.5);
  L_corr  ~ lkj_corr_cholesky(2);
  for (m in 1:M)
    b[m, ] ~ multi_normal_cholesky(rep_vector(0, 2),
               diag_pre_multiply(sigma_b, L_corr));
  y ~ normal(mu, sigma_e);
}

generated quantities {
  matrix[2, 2] Omega;
  Omega = multiply_lower_tri_self_transpose(L_corr);
  real rho = Omega[1, 2];
  real icc = sigma_b[1]^2 / (sigma_b[1]^2 + sigma_e^2);
  vector[N] log_lik;
  vector[N] y_rep;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | mu[n], sigma_e);
    y_rep[n]   = normal_rng(mu[n], sigma_e);
  }
}
"

writeLines(stan_model_code, "aids_lmm.stan")

# ---------------------------------------------------------
# 5. Fit model
# ---------------------------------------------------------

set.seed(123)

fit <- stan(
  file    = "aids_lmm.stan",
  data    = stan_data,
  chains  = 4,
  iter    = 6000,
  warmup  = 1000,
  seed    = 123,
  control = list(adapt_delta   = 0.95,
                 max_treedepth = 12)
)

# ---------------------------------------------------------
# 6. Posterior summaries
# ---------------------------------------------------------

print(fit,
      pars   = c("beta1","beta2","beta3",
                 "beta4","beta5","beta6",
                 "sigma_b","sigma_e","rho","icc"),
      digits = 3)

# ---------------------------------------------------------
# 7. Clean results table
# ---------------------------------------------------------

post <- extract(fit)

param_list <- list(
  beta1    = post$beta1,
  beta2    = post$beta2,
  beta3    = post$beta3,
  beta4    = post$beta4,
  beta5    = post$beta5,
  beta6    = post$beta6,
  sigma_b0 = post$sigma_b[, 1],
  sigma_b1 = post$sigma_b[, 2],
  sigma_e  = post$sigma_e,
  rho      = post$rho,
  icc      = post$icc
)

param_labels <- c(
  "$\\beta_1$ (intercept)",
  "$\\beta_2$ (time)",
  "$\\beta_3$ (drug $\\times$ time)",
  "$\\beta_4$ (gender)",
  "$\\beta_5$ (prior OI)",
  "$\\beta_6$ (AZT)",
  "$\\sigma_{b_0}$ (intercept SD)",
  "$\\sigma_{b_1}$ (slope SD)",
  "$\\sigma_e$ (residual SD)",
  "$\\rho$ (intercept--slope corr.)",
  "ICC"
)

results <- data.frame(
  Parameter = param_labels,
  Mean      = sapply(param_list, mean),
  SD        = sapply(param_list, sd),
  Lower95   = sapply(param_list, quantile, 0.025),
  Upper95   = sapply(param_list, quantile, 0.975),
  row.names = NULL
)
results[, 2:5] <- round(results[, 2:5], 3)
print(results, row.names = FALSE)

# ---------------------------------------------------------
# 8. Trace plots
# ---------------------------------------------------------

chain_cols <- c("steelblue","tomato","forestgreen","purple")

sims     <- extract(fit, permuted = FALSE)
n_iter   <- dim(sims)[1]
n_chains <- dim(sims)[2]

to_mcmc_list <- function(param) {
  as.mcmc.list(
    lapply(seq_len(n_chains),
           function(j) mcmc(sims[, j, param]))
  )
}

library(coda)

beta1_mcmc   <- to_mcmc_list("beta1")
beta2_mcmc   <- to_mcmc_list("beta2")
beta5_mcmc   <- to_mcmc_list("beta5")
sigmab_mcmc  <- to_mcmc_list("sigma_b[1]")
sigmae_mcmc  <- to_mcmc_list("sigma_e")
rho_mcmc     <- to_mcmc_list("rho")

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
# 9. Autocorrelation plots
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
# 10. Gelman-Rubin diagnostics
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
# 11. Effective sample sizes
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
# 12. Posterior density plots
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
# 13. Posterior predictive check
# ---------------------------------------------------------

y_rep <- post$y_rep

p_ppc <- ppc_dens_overlay(df$sqrtCD4, y_rep[1:100, ]) +
  labs(title = "Posterior predictive check",
       x     = expression(sqrt(CD4))) +
  theme_bw(base_size = 11)
print(p_ppc)
ggsave("Plots/aids_ppc.png", p_ppc,
       width = 7, height = 4, dpi = 300)

# ---------------------------------------------------------
# 14. LOO cross-validation
# ---------------------------------------------------------

log_lik    <- extract_log_lik(fit, merge_chains = FALSE)
r_eff      <- relative_eff(exp(log_lik))
loo_result <- loo(log_lik, r_eff = r_eff)
print(loo_result)