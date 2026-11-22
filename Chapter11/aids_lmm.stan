
data {
  int<lower=1> N;
  int<lower=1> M;
  vector[N] y;
  vector[N] time;
  array[N] int<lower=0, upper=1> drug;
  array[N] int<lower=0, upper=1> gender;
  array[N] int<lower=0, upper=1> prevOI;
  array[N] int<lower=0, upper=1> AZT;
  array[N] int<lower=1, upper=M> patient;
}

parameters {
  real beta1;
  real beta2;
  real beta3;
  real beta4;
  real beta5;
  real beta6;
  matrix[M, 2] b;
  cholesky_factor_corr[2] L_corr;
  vector<lower=0>[2] sigma_b;
  real<lower=0> sigma_e;
}

transformed parameters {
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
  beta1 ~ normal(0, 10);
  beta2 ~ normal(0, 10);
  beta3 ~ normal(0, 10);
  beta4 ~ normal(0, 10);
  beta5 ~ normal(0, 10);
  beta6 ~ normal(0, 10);
  sigma_b ~ cauchy(0, 2.5);
  sigma_e ~ cauchy(0, 2.5);
  L_corr  ~ lkj_corr_cholesky(2);
  for (m in 1:M)
    b[m, ] ~ multi_normal_cholesky(rep_vector(0, 2),
               diag_pre_multiply(sigma_b, L_corr));
  y ~ normal(mu, sigma_e);
}

generated quantities {
  matrix[2, 2] Omega;
  Omega = multiply_lower_tri_self_transpose(L_corr);
  real rho = Omega[1, 2];
  real icc = sigma_b[1]^2 / (sigma_b[1]^2 + sigma_e^2);
  vector[N] log_lik;
  vector[N] y_rep;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | mu[n], sigma_e);
    y_rep[n]   = normal_rng(mu[n], sigma_e);
  }
}

