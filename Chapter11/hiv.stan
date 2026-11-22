data {
  int<lower=1> N;                          // total observations
  int<lower=1> M;                          // number of patients
  vector[N] y;                             // sqrt(CD4)
  vector[N] time;                          // observation time (months)
  array[N] int<lower=0, upper=1> drug;     // 1 = ddI, 0 = ddC
  array[N] int<lower=0, upper=1> gender;   // gender indicator
  array[N] int<lower=0, upper=1> prevOI;   // prior opportunistic infection
  array[N] int<lower=0, upper=1> AZT;      // prior AZT use
  array[N] int<lower=1, upper=M> patient;  // patient index
}

parameters {
  // fixed effects
  real beta1;                              // intercept
  real beta2;                              // time effect
  real beta3;                              // drug x time interaction
  real beta4;                              // gender
  real beta5;                              // prior OI
  real beta6;                              // AZT stratum
  
  // random effects
  matrix[M, 2] b;                          // [b0i, b1i] per patient
  
  // random effects covariance (Cholesky factor)
  cholesky_factor_corr[2] L_corr;         // Cholesky of correlation matrix
  vector<lower=0>[2] sigma_b;             // random effect SDs
  
  // residual SD
  real<lower=0> sigma_e;
}

transformed parameters {
  // reconstruct covariance matrix
  matrix[2, 2] D;
  D = diag_pre_multiply(sigma_b, L_corr) *
    diag_pre_multiply(sigma_b, L_corr)';

  vector[N] mu;
  for (n in 1:N) {
    mu[n] = beta1
          + beta2  * time[n]
          + beta3  * (time[n] * drug[n])
          + beta4  * gender[n]
          + beta5  * prevOI[n]
          + beta6  * AZT[n]
          + b[patient[n], 1]
          + b[patient[n], 2] * time[n];
  }
}

model {
  // priors on fixed effects
  beta1 ~ normal(0, 10);
  beta2 ~ normal(0, 10);
  beta3 ~ normal(0, 10);
  beta4 ~ normal(0, 10);
  beta5 ~ normal(0, 10);
  beta6 ~ normal(0, 10);

  // priors on variance components
  sigma_b ~ cauchy(0, 2.5);
  sigma_e ~ cauchy(0, 2.5);

  // LKJ prior on correlation matrix
  L_corr ~ lkj_corr_cholesky(2);

  // random effects — non-centred parameterisation
  for (m in 1:M) {
    b[m, ] ~ multi_normal_cholesky(
               rep_vector(0, 2),
               diag_pre_multiply(sigma_b, L_corr)
             );
  }

  // likelihood
  y ~ normal(mu, sigma_e);
}

generated quantities {
  // correlation between random intercept and slope
  matrix[2, 2] Omega;
  Omega = multiply_lower_tri_self_transpose(L_corr);
  real rho = Omega[1, 2];

  // ICC at time 0
  real icc = sigma_b[1]^2 /
             (sigma_b[1]^2 + sigma_e^2);

  // pointwise log-likelihood for LOO
  vector[N] log_lik;
  vector[N] y_rep;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | mu[n], sigma_e);
    y_rep[n]   = normal_rng(mu[n], sigma_e);
  }
}