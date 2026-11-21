data {
  int<lower=1> N;
  vector<lower=0>[N] y;    // income.1
  int<lower=1> K;
  matrix[N, K] X;
}

parameters {
  vector[K] beta;
  real<lower=0> sigma;     // standard deviation of log(y)
}

model {
  vector[N] mu = X * beta; // meanlog

  beta ~ normal(0,100);
  sigma ~ exponential(1);

  y ~ lognormal(mu, sigma);
}
