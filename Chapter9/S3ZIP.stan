data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> Q;
  int<lower=0> Y[N];
  matrix[N,K] X;
  matrix[N,Q] Z;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_gamma;
}

parameters {
  vector[K] beta;
  vector[Q] gamma;
}

transformed parameters {
  vector[N] lambda = exp(X * beta);
  vector[N] kappa = inv_logit(Z * gamma);
}

model {
  beta ~ normal(0, sigma_beta);
  gamma ~ normal(0, sigma_gamma);

  for (n in 1:N) {
    if (Y[n] == 0) {
      target += log_mix(kappa[n], 0, poisson_lpmf(0 | lambda[n]));
    } else {
      target += log1m(kappa[n]) + poisson_lpmf(Y[n] | lambda[n]);
    }
  }
}

generated quantities {
  vector[N] log_lik;
  int y_rep[N];

  for (n in 1:N) {
    if (Y[n] == 0) {
      log_lik[n] = log_mix(kappa[n], 0, poisson_lpmf(0 | lambda[n]));
    } else {
      log_lik[n] = log1m(kappa[n]) + poisson_lpmf(Y[n] | lambda[n]);
    }

    if (bernoulli_rng(kappa[n]) == 1) {
      y_rep[n] = 0;
    } else {
      y_rep[n] = poisson_rng(lambda[n]);
    }
  }
}
