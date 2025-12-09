data {
  int<lower=1> N;               // number of observations
  vector<lower=0>[N] y;         // response variable (income)
  int<lower=1> K;               // number of covariates
  matrix[N, K] X;               // design matrix
}

parameters {
  vector[K] beta;               // regression coefficients
  real<lower=0> alpha;          // shape parameter
}

transformed parameters {
  vector[N] mu;                 // mean of Gamma
  mu = exp(X * beta);           // log link: mu_i = exp(X_i * beta)
}

model {
  // Priors
  beta ~ normal(0, 5);
  alpha ~ exponential(1);

  // Likelihood
  for (n in 1:N)
    y[n] ~ gamma(alpha, mu[n] / alpha);  // scale = mu / alpha
}
