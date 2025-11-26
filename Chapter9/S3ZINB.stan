data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> Q;
  int<lower=0> Y[N];
  matrix[N,K] X;
  matrix[N,Q] Z;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_gamma;
  real<lower=0> a_phi;
  real<lower=0> b_phi;
}

parameters {
  vector[K] beta;
  vector[Q] gamma;
  real<lower=0> phi;
}

transformed parameters {
  vector[N] lambda = exp(X * beta);
  vector[N] kappa = inv_logit(Z * gamma);
}

model {
  beta ~ normal(0, sigma_beta);
  gamma ~ normal(0, sigma_gamma);
  phi ~ gamma(a_phi, b_phi);

  for (n in 1:N) {
    if (Y[n] == 0) {
      target += log_mix(kappa[n], 0, neg_binomial_2_lpmf(0 | lambda[n], phi));
    } else {
      target += log1m(kappa[n]) + neg_binomial_2_lpmf(Y[n] | lambda[n], phi);
    }
  }
}

generated quantities {
  vector[N] log_lik;
  int y_rep[N];

  for (n in 1:N) {
    if (Y[n] == 0) {
      log_lik[n] = log_mix(kappa[n], 0, neg_binomial_2_lpmf(0 | lambda[n], phi));
    } else {
      log_lik[n] = log1m(kappa[n]) + neg_binomial_2_lpmf(Y[n] | lambda[n], phi);
    }

    if (bernoulli_rng(kappa[n]) == 1) {
      y_rep[n] = 0;
    } else {
      y_rep[n] = neg_binomial_2_rng(lambda[n], phi);
    }
  }
}
