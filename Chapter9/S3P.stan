data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=0> Y[N];
  matrix[N,K] X;
}

parameters {
  vector[K] beta;
}

transformed parameters {
  vector[N] lambda = exp(X * beta);
}

model {
  beta ~ normal(0, 100);
  Y ~ poisson(lambda);
}

generated quantities {
  vector[N] log_lik;

  for (n in 1:N) {
    log_lik[n] = poisson_lpmf(Y[n] | lambda[n]);
  }
}
