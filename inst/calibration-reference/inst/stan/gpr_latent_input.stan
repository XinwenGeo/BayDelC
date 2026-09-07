// Latent-input Gaussian Process sensitivity model for BayDelC.
// x_obs is standardized O2. x_error is standardized O2 uncertainty.
// y_obs is Δδ13C in per mil. y_error is response-side uncertainty.

data {
  int<lower=1> N;
  vector[N] x_obs;
  vector<lower=0>[N] x_error;
  vector[N] y_obs;
  vector<lower=0>[N] y_error;
  real y_bar;
}
parameters {
  vector[N] x_true;
  real alpha;
  real<lower=0> eta;
  real<lower=0> rho;
  real<lower=0> sigma_model;
  vector[N] z;
}
transformed parameters {
  matrix[N, N] K;
  matrix[N, N] L_K;
  vector[N] f_latent;
  vector[N] mu;

  for (i in 1:N) {
    for (j in i:N) {
      K[i, j] = square(eta) * exp(-0.5 * square((x_true[i] - x_true[j]) / rho));
      K[j, i] = K[i, j];
    }
  }
  for (i in 1:N) K[i, i] = K[i, i] + 1e-8;
  L_K = cholesky_decompose(K);
  f_latent = L_K * z;
  mu = alpha + f_latent;
}
model {
  alpha ~ normal(y_bar, 2);
  eta ~ exponential(4);                         // shrinkage on GP amplitude
  rho ~ lognormal(log(0.8), 0.4);               // smoothness on standardized O2 scale
  sigma_model ~ exponential(2);
  z ~ normal(0, 1);
  x_true ~ normal(x_obs, x_error);

  for (n in 1:N) {
    y_obs[n] ~ normal(mu[n], sqrt(square(y_error[n]) + square(sigma_model)));
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  for (n in 1:N) {
    real sd_n = sqrt(square(y_error[n]) + square(sigma_model));
    log_lik[n] = normal_lpdf(y_obs[n] | mu[n], sd_n);
    y_rep[n] = normal_rng(mu[n], sd_n);
  }
}
