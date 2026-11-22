data {
  int<lower=1> N;                        // total number of observations
  int<lower=1> M;                        // number of cows
  int<lower=1> K;                        // number of fixed effects
  vector[N] y;                           // milk protein content
  matrix[N, K] X;                        // fixed-effects design matrix
  array[N] int<lower=1, upper=M> cow;    // cow index for each observation
  vector[N] time;                        // week of measurement
}

parameters {
  vector[K] beta;                        // fixed effects
  vector[M] b;                           // cow-level random intercepts
  real<lower=0> tau_b;
  real<lower=0> tau_e;
}

transformed parameters {
  real<lower=0> sigma_b = 1.0 / sqrt(tau_b);
  real<lower=0> sigma_e = 1.0 / sqrt(tau_e);
}


model {
  // priors
  beta    ~ normal(0, 2.5);
  tau_b ~ gamma(0.001, 0.001);
  tau_e ~ gamma(0.001, 0.001);

  // random effects
  b ~ normal(0, sigma_b);

  // likelihood
  y ~ normal(X * beta + b[cow], sigma_e);
}

generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  real icc;

  icc = sigma_b^2 / (sigma_b^2 + sigma_e^2);

  for (n in 1:N) {
    real mu_n  = X[n] * beta + b[cow[n]];
    log_lik[n] = normal_lpdf(y[n] | mu_n, sigma_e);
    y_rep[n]   = normal_rng(mu_n, sigma_e);
  }
}