data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=0> Y[N];
  matrix[N,K] X;
}

parameters {
  vector[K] beta;
  real<lower=0> phi;  
}

transformed parameters {
  vector[N] lambda = exp(X * beta);
}

model {
  beta ~ normal(0, 100);
  phi ~ gamma(0.1, 0.1);
  Y ~ neg_binomial_2(lambda, phi);
}




generated quantities {
  vector[N] log_lik;

  for (n in 1:N) {
    log_lik[n] = neg_binomial_2_lpmf(Y[n] | lambda[n], phi);
  }
}
