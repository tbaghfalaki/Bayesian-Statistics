data {
  int<lower=1> N;                  // number of observations
  vector[N] logT;                  // log survival times
  int<lower=0,upper=1> delta[N];   // event indicator (1=event, 0=censored)
  int<lower=1> K;                  // number of covariates
  matrix[N, K] X;                  // design matrix
}

parameters {
  vector[K] beta;
  real<lower=0> sigma;             // scale parameter
}

model {
  // Priors
  beta ~ normal(0, 100);
  sigma ~ cauchy(0, 2.5);

  // Likelihood
  for (i in 1:N) {
    if (delta[i] == 1) {
      target += normal_lpdf(logT[i] | X[i] * beta, sigma);
    } else {
      target += normal_lccdf(logT[i] | X[i] * beta, sigma);
    }
  }
}

generated quantities {
  vector[N] log_lik;

  for (i in 1:N) {
    if (delta[i] == 1) {
      log_lik[i] = normal_lpdf(logT[i] | dot_product(X[i], beta), sigma);
    } else {
      log_lik[i] = normal_lccdf(logT[i] | dot_product(X[i], beta), sigma);
    }
  }
}

