# =========================================================
# Bayesian Logistic Regression: Survival ~ Predictors
# =========================================================

rm(list = ls())

library(R2jags)
library(coda)

# ---------------------------------------------------------
# 1. Data preparation
# ---------------------------------------------------------

Data <- read.csv(
  "/Users/user/Downloads/Bayesian-Statistics-main/Chapter1/train_complete.csv",
  stringsAsFactors = FALSE
)

cat("Dimensions:", dim(Data), "\n")
cat("Variables:", names(Data), "\n")
head(Data)

# ---------------------------------------------------------
# Data cleaning & validation
# ---------------------------------------------------------

# Inspect Survived
cat("Survived values (unique):", unique(Data$Survived), "\n")
cat("Survived NAs:", sum(is.na(Data$Survived)), "\n")

# Force Survived to integer 0/1
Data$Survived <- as.integer(Data$Survived)

# Check: must be exactly 0 or 1
bad <- !Data$Survived %in% c(0L, 1L)
if (any(bad, na.rm = TRUE))
  cat("WARNING: unexpected Survived values at rows:", which(bad), "\n")

# Drop rows where Survived is NA or not 0/1,
# or where any predictor is NA
keep <- Data$Survived %in% c(0L, 1L) &
  !is.na(Data$Sex)      &
  !is.na(Data$Age)      &
  !is.na(Data$Pclass)   &
  !is.na(Data$Fare)     &
  !is.na(Data$Embarked) &
  Data$Embarked != ""   &   # drop blank Embarked levels
  !is.na(Data$SibSp)    &
  !is.na(Data$Parch)

cat("Rows kept:", sum(keep), "| Rows dropped:", sum(!keep), "\n")
Data <- Data[keep, ]

# Scale continuous predictors to avoid extreme inprod values at init
Data$Age   <- scale(Data$Age)
Data$Fare  <- scale(Data$Fare)
Data$SibSp <- scale(Data$SibSp)
Data$Parch <- scale(Data$Parch)

# Design matrix (intercept included automatically)
X <- model.matrix(
  ~ Sex + Age + as.factor(Pclass) + Fare +
    as.factor(Embarked) + SibSp + Parch,
  data = Data
)

K <- ncol(X)
cat("Number of predictors (incl. intercept):", K, "\n")
cat("Column names:\n"); print(colnames(X))

# Sanity check: no NA or Inf in X
stopifnot(!any(is.na(X)), !any(is.infinite(X)))

data_jags <- list(
  N        = nrow(X),
  K        = K,
  Survived = Data$Survived,
  X        = X
)

# ---------------------------------------------------------
# 2. JAGS model
# ---------------------------------------------------------

logistic_model <- "
model {

  # Likelihood
  for (i in 1:N) {
    Survived[i] ~ dbern(pi[i])
    logit(pi[i]) <- inprod(beta[], X[i, ])
  }

  # Weakly informative priors
  for (k in 1:K) {
    beta[k] ~ dnorm(0, 0.001)
  }
}
"

# ---------------------------------------------------------
# 3. MCMC settings
# ---------------------------------------------------------

# Diagnose row 24
cat("\nRow 24 of X:\n");  print(X[24, ])
cat("Survived[24]:", Data$Survived[24], "\n")

# Very small inits so inprod(beta, X[i,]) stays near 0 -> pi near 0.5
inits_fun <- function() {
  list(beta = rep(0, K))
}

params <- "beta"


# ---------------------------------------------------------
# 4. Model fitting
# ---------------------------------------------------------

set.seed(42)

fit <- jags(
  data               = data_jags,
  inits              = inits_fun,
  parameters.to.save = params,
  model.file         = textConnection(logistic_model),
  n.chains           = 3,
  n.iter             = 15000,
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

sims     <- fit$BUGSoutput$sims.array   # [iterations x chains x parameters]
n_iter   <- dim(sims)[1]
n_chains <- dim(sims)[2]

# Beta parameter names from design matrix
beta_names <- colnames(X)   # e.g. "(Intercept)", "Sexmale", "Age", ...
beta_jags  <- paste0("beta[", seq_len(K), "]")   # names as stored in sims

# Helper: build mcmc.list for one parameter
to_mcmc_list <- function(param) {
  as.mcmc.list(
    lapply(seq_len(n_chains), function(j) mcmc(sims[, j, param]))
  )
}

all_mcmc <- setNames(
  lapply(beta_jags, to_mcmc_list),
  beta_names
)

# ---------------------------------------------------------
# 7. Trace plots  — 3x3 panels, one figure per 9 parameters
# ---------------------------------------------------------

chain_cols <- c("steelblue", "tomato", "forestgreen")

trace_param <- function(mcmc_list, label) {
  plot(
    NA,
    xlim = c(1, n_iter),
    ylim = range(unlist(mcmc_list)),
    type = "n",
    main = label,
    xlab = "Iteration",
    ylab = ""
  )
  for (j in seq_len(n_chains))
    lines(as.numeric(mcmc_list[[j]]), col = chain_cols[j], lwd = 0.6)
}

param_names <- names(all_mcmc)
pages_idx   <- split(seq_along(param_names),
                     ceiling(seq_along(param_names) / 9))

for (pg in pages_idx) {
  par(mfrow = c(3, 3), mar = c(3, 3, 2, 1), oma = c(0, 0, 3, 0))
  for (i in pg) trace_param(all_mcmc[[i]], param_names[i])
  mtext("Trace plots (all chains)", outer = TRUE, cex = 1.2, font = 2)
  legend("topright", legend = paste("Chain", seq_len(n_chains)),
         col = chain_cols, lty = 1, bty = "n", cex = 0.9, xpd = NA)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0))
}

# ---------------------------------------------------------
# 8. Autocorrelation plots — chain 1, 3x3 panels
# ---------------------------------------------------------

acf_vals <- lapply(all_mcmc, function(m) autocorr(m[[1]]))
lags      <- as.numeric(dimnames(acf_vals[[1]])[[1]])

for (pg in pages_idx) {
  par(mfrow = c(3, 3), mar = c(4, 4, 2, 1), oma = c(0, 0, 3, 0))
  for (i in pg) {
    av <- acf_vals[[i]][, 1, 1]
    barplot(av, names.arg = lags, main = param_names[i],
            xlab = "Lag", ylab = "ACF",
            col = "steelblue", border = NA, ylim = c(-1, 1))
    abline(h = 0,             lty = 1)
    abline(h = c(-0.1, 0.1), lty = 2, col = "tomato")
  }
  mtext("Autocorrelation plots - Chain 1", outer = TRUE, cex = 1.2, font = 2)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0))
}

# ---------------------------------------------------------
# 9. Gelman-Rubin diagnostic  (R-hat < 1.1 is good)
# ---------------------------------------------------------

cat("\n--- Gelman-Rubin diagnostics ---\n")
gr <- lapply(all_mcmc, gelman.diag)
lapply(seq_along(gr), function(i) {
  cat(param_names[i], ":\n"); print(gr[[i]])
})

# ---------------------------------------------------------
# 10. Effective sample size
# ---------------------------------------------------------

cat("\n--- Effective sample sizes ---\n")
ess <- sapply(all_mcmc, effectiveSize)
names(ess) <- param_names
print(round(ess))

# ---------------------------------------------------------
# 11. Posterior density plots — 3x3 panels
# ---------------------------------------------------------

pool <- function(m) unlist(lapply(m, as.numeric))
post <- lapply(all_mcmc, pool)

for (pg in pages_idx) {
  par(mfrow = c(3, 3), mar = c(4, 4, 2, 1), oma = c(0, 0, 3, 0))
  for (i in pg) {
    nm <- param_names[i]
    d  <- density(post[[nm]])
    ci <- quantile(post[[nm]], c(0.025, 0.975))
    
    plot(d, main = nm, xlab = "", lwd = 2, col = "steelblue")
    
    x_shade <- c(ci[1], d$x[d$x >= ci[1] & d$x <= ci[2]], ci[2])
    y_shade <- c(0,     d$y[d$x >= ci[1] & d$x <= ci[2]], 0)
    polygon(x_shade, y_shade, col = adjustcolor("steelblue", 0.25), border = NA)
    
    abline(v = mean(post[[nm]]), lty = 2, col = "tomato", lwd = 1.5)
    legend("topright",
           legend = c(sprintf("Mean = %.3f", mean(post[[nm]])),
                      sprintf("95%% CI [%.2f, %.2f]", ci[1], ci[2])),
           bty = "n", cex = 0.65)
  }
  mtext("Posterior density plots (95% CI shaded)", outer = TRUE, cex = 1.2, font = 2)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0))
}

# ---------------------------------------------------------
# 12. Odds ratios (exp(beta)) with 95% CI
# ---------------------------------------------------------

cat("\n--- Odds Ratios (posterior mean + 95% CI) ---\n")
or_summary <- t(sapply(post, function(p) {
  ep <- exp(p)
  c(OR    = mean(ep),
    Lower = quantile(ep, 0.025),
    Upper = quantile(ep, 0.975))
}))
rownames(or_summary) <- param_names
print(round(or_summary, 3))

