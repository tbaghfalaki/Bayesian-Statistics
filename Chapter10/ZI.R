rm(list = ls())
library(bellreg)
library(rstan)
setwd('/Users/user/Desktop/Lab/Example3')  # change if necessary
# 1. Data -------------------------------------------------------------------
data(cells, package = 'bellreg')
cells$weightf <- factor(as.character(cells$weight),
                        levels = c('normal', 'over', 'obese'))
cells$agef <- factor(as.character(cells$age),
                     levels = c('young', 'mid', 'old'))
cells$smoker <- factor(cells$smoker)
cells$gender <- factor(cells$gender)
Y <- as.integer(cells$cells)
X <- model.matrix(~ smoker + weightf + agef + gender, data = cells)
Z <- X
N <- length(Y); K <- ncol(X); Q <- ncol(Z)
cat('N =', N, '; K =', K, '; Q =', Q, '\n')
print(colnames(X))

barplot(table(factor(Y, levels = 0:max(Y))), col = 'skyblue',
        border = 'white', main = 'Observed infected-cell counts',
        xlab = 'Number of infected cells', ylab = 'Frequency')

# 2. Fit models -------------------------------------------------------------
rstan_options(auto_write = TRUE)
stan_data <- list(N = N, K = K, Q = Q, Y = Y, X = X, Z = Z,
                  sigma_beta = 2.5, sigma_gamma = 2.5)
zinb_data <- c(stan_data, list(a_phi = 0.1, b_phi = 0.1))

ZIP  <- rstan::stan_model('S3ZIP.stan')
ZINB <- rstan::stan_model('S3ZINB.stan')
fit_zip  <- rstan::sampling(ZIP,  data = stan_data, chains = 2, cores = 2,
                            iter = 3000, warmup = 1000, seed = 123,
                            control = list(adapt_delta = 0.95))
fit_zinb <- rstan::sampling(ZINB, data = zinb_data, chains = 2, cores = 2,
                            iter = 3000, warmup = 1000, seed = 123,
                            control = list(adapt_delta = 0.95))
print(fit_zip,  pars = c('beta', 'gamma'))
print(fit_zinb, pars = c('beta', 'gamma', 'phi'))

# 3. Posterior draws --------------------------------------------------------
zip  <- rstan::extract(fit_zip,  pars = c('beta','gamma','log_lik','y_rep'))
zinb <- rstan::extract(fit_zinb, pars = c('beta','gamma','phi','log_lik','y_rep'))

# 4. Criteria ---------------------------------------------------------------
ll_zip  <- mean(rowSums(zip$log_lik))
ll_zinb <- mean(rowSums(zinb$log_lik))

BIC_zip  <- -2 * ll_zip  + (K + Q)      * log(N)
BIC_zinb <- -2 * ll_zinb + (K + Q + 1L) * log(N)

logmeanexp <- function(x) { m <- max(x); m + log(mean(exp(x - m))) }
LPML_zip  <- sum(apply(zip$log_lik,  2, function(x) -logmeanexp(-x)))
LPML_zinb <- sum(apply(zinb$log_lik, 2, function(x) -logmeanexp(-x)))

DIC_zip  <- -4 * mean(rowSums(zip$log_lik))  + 2 * ll_zip
DIC_zinb <- -4 * mean(rowSums(zinb$log_lik)) + 2 * ll_zinb

criteria <- data.frame(
  Model = c('ZIP', 'ZINB'),
  BIC   = c(BIC_zip,   BIC_zinb),
  LPML  = c(LPML_zip,  LPML_zinb),
  DIC   = c(DIC_zip,   DIC_zinb)
)
print(criteria, digits = 8, row.names = FALSE)

