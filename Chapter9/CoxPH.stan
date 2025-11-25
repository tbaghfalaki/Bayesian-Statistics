data {
  int<lower=1> N;                  // number of observations
  int<lower=1> K;                  // number of covariates
  matrix[N, K] X;                  // design matrix (no intercept)
  vector<lower=0>[N] T;            // survival times
  int<lower=0,upper=1> delta[N];   // event indicator
}

parameters {
  vector[K] beta;                  // regression coefficients
  real<lower=0> lambda0;           // constant baseline hazard
}

model {
  // Priors
  beta ~ normal(0, 5);
  lambda0 ~ gamma(1, 1);

  // Likelihood
  for (i in 1:N) {
    real hazard = lambda0 * exp(dot_product(X[i], beta));
    if (delta[i] == 1) {
      target += log(hazard) - hazard * T[i];
    } else {
      target += -hazard * T[i];
    }
  }
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    real hazard = lambda0 * exp(dot_product(X[i], beta));
    if (delta[i] == 1) {
      log_lik[i] = log(hazard) - hazard * T[i];
    } else {
      log_lik[i] = -hazard * T[i];
    }
  }
}
