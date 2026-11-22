data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> K;
  matrix[N, K] X;
  array[N] int<lower=1, upper=J> school;
  vector[N] y;
}

parameters {
  vector[K] beta;
  vector[J] b;
  real<lower=0> tau_b;
  real<lower=0> tau_e;
}

transformed parameters {
  real<lower=0> sigma_b = 1.0 / sqrt(tau_b);
  real<lower=0> sigma_e = 1.0 / sqrt(tau_e);
}

model {
  beta  ~ normal(0, 1000);
  tau_b ~ gamma(0.001, 0.001);
  tau_e ~ gamma(0.001, 0.001);

  b ~ normal(0, sigma_b);

  y ~ normal(X * beta + b[school], sigma_e);
}

generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  real icc;

  icc = sigma_b^2 / (sigma_b^2 + sigma_e^2);

  for (n in 1:N) {
    real mu_n  = X[n] * beta + b[school[n]];
    log_lik[n] = normal_lpdf(y[n] | mu_n, sigma_e);
    y_rep[n]   = normal_rng(mu_n, sigma_e);
  }
}

