data {
  int<lower=1> N;              // number of observations
  int<lower=2> K;              // number of ordinal categories
  int<lower=1> P;              // number of covariates
  matrix[N, P] X;              // design matrix
  int<lower=1, upper=K> y[N]; // ordinal outcome (1..K)
}

parameters {
  vector[P] beta;              // regression coefficients
  ordered[K-1] c;              // thresholds (cutpoints)
}

model {
  // Priors
  beta ~ normal(0, 10);         // weakly informative prior
  c ~ normal(0, 10);            // prior for cutpoints

  // Likelihood
  for (n in 1:N)
    y[n] ~ ordered_logistic(X[n] * beta, c);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = ordered_logistic_lpmf(y[n] | X[n] * beta, c);
}
