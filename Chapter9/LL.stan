data {
  int<lower=1> N;                  // number of observations
  vector<lower=0>[N] T;         // survival times  
  int<lower=0,upper=1> delta[N];   // event indicator (1=event, 0=censored)
  int<lower=1> K;                  // number of covariates
  matrix[N, K] X;                  // design matrix
}

parameters {
  vector[K] beta;
  real<lower=0> sigma;             // scale parameter of logistic error
}

model {
  vector[N] eta = X * beta;

  // Priors
  beta ~ normal(0, 100);
  sigma ~ cauchy(0, 2.5);

  // Likelihood
  for (i in 1:N) {
    real z = (log(T[i]) - eta[i]) / sigma;

    if (delta[i] == 1) {
      target += -log(sigma) - log(T[i]) + logistic_lpdf(z | 0, 1);
    } else {
      target += logistic_lccdf(z | 0, 1);
    }
  }
}

generated quantities {
  vector[N] log_lik;
  
  for (i in 1:N) {
    real z = (log(T[i]) - dot_product(X[i], beta)) / sigma;

    if (delta[i] == 1) {
      log_lik[i] = -log(sigma) - log(T[i]) + logistic_lpdf(z | 0, 1);
    } else {
      log_lik[i] = logistic_lccdf(z | 0, 1);
    }
  }
}

