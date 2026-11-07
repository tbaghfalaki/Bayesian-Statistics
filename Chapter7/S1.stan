// model1.stan
data {
  int<lower=0> N;         // number of observations
  vector[N] x;            // predictor (sex, coded as 0/1)
  vector[N] y;            // outcome (BMI)
}

parameters {
  real beta0;             // intercept
  real beta1;             // slope
  real<lower=0> sigma;    // standard deviation
}

model {
  // Priors
  beta0 ~ normal(0, 1e3);
  beta1 ~ normal(0, 1e3);
  sigma ~ cauchy(0, 5);  // weakly informative prior for SD
  
  // Likelihood
  y ~ normal(beta0 + beta1 * x, sigma);
}

generated quantities {
  vector[N] mu;
  for (i in 1:N)
    mu[i] = beta0 + beta1 * x[i];
}
