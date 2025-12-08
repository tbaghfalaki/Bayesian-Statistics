data {
  int<lower=1> N;                  
  vector<lower=0>[N] T;            
  int<lower=0,upper=1> delta[N];   
  int<lower=1> K;                  
  matrix[N, K] X;                  
}

parameters {
  vector[K] beta;                  
  real<lower=0> alpha;             
}

model {
  beta ~ normal(0, 100);
  alpha ~ gamma(1, 1);

  for (i in 1:N) {
    real lambda_i = exp(-X[i] * beta);   // lambda absorbed in beta
    if (delta[i] == 1)
      target += weibull_lpdf(T[i] | alpha, lambda_i);
    else
      target += weibull_lccdf(T[i] | alpha, lambda_i);
  }
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    real lambda_i = exp(-X[i] * beta);
    if (delta[i] == 1)
      log_lik[i] = weibull_lpdf(T[i] | alpha, lambda_i);
    else
      log_lik[i] = weibull_lccdf(T[i] | alpha, lambda_i);
  }
}

