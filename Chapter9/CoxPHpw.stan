data {
  int<lower=1> N;                  // number of subjects
  int<lower=1> K;                  // number of covariates
  int<lower=1> J;                  // number of intervals
  matrix[N,K] X;                    // design matrix
  vector<lower=0>[N] t;             // survival times
  int<lower=0,upper=1> delta[N];    // censoring indicators
  int<lower=1,upper=J> interval[N]; // interval index for each subject
  vector[J+1] breaks;               // break points for intervals
}

parameters {
  vector[K] beta;                   // regression coefficients
  vector<lower=0>[J] lambda;        // piecewise baseline hazards
}

model {
  // Priors
  beta ~ normal(0,5);
  lambda ~ gamma(1,1);

  // Likelihood
  for (i in 1:N) {
    real hazard = lambda[interval[i]] * exp(dot_product(X[i], beta));
    real cum_hazard = 0;
    
    // Sum of hazards over completed intervals
    for (j in 1:(interval[i]-1))
      cum_hazard += lambda[j] * (breaks[j+1] - breaks[j]) * exp(dot_product(X[i], beta));
    
    // Partial interval
    cum_hazard += lambda[interval[i]] * (t[i] - breaks[interval[i]]) * exp(dot_product(X[i], beta));

    if (delta[i] == 1)
      target += log(hazard) - cum_hazard;
    else
      target += -cum_hazard;
  }
}



generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    real hazard = lambda[interval[i]] * exp(dot_product(X[i], beta));
    real cum_hazard = 0;

    for (j in 1:(interval[i]-1))
      cum_hazard += lambda[j] * (breaks[j+1] - breaks[j]) * exp(dot_product(X[i], beta));

    cum_hazard += lambda[interval[i]] * (t[i] - breaks[interval[i]]) * exp(dot_product(X[i], beta));

    if (delta[i] == 1)
      log_lik[i] = log(hazard) - cum_hazard;
    else
      log_lik[i] = -cum_hazard;
  }
}
