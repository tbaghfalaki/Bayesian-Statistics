data {
  int<lower=1> N;
  vector<lower=0>[N] y;    // income.1
  int<lower=1> K;
  matrix[N, K] X;
}

parameters {
  vector[K] beta;
}

model {
  vector[N] eta = X * beta;
  beta ~ normal(0, 100);
  y ~ exponential(exp(-eta));   // rate = exp(eta)
}
