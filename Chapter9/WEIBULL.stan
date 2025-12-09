data {
  int<lower=1> N;
  vector<lower=0>[N] y;    // income.1
  int<lower=1> K;
  matrix[N, K] X;
}

parameters {
  vector[K] beta;
  real<lower=0> shape;   // k parameter
}

model {
  vector[N] scale = exp(X * beta);  // scale parameter

  beta ~ normal(0, 5);
  shape ~ normal(2, 1);

  for (n in 1:N)
    y[n] ~ weibull(shape, scale[n]);
}
