data {
  int<lower=1> N;
  vector[N] x_obs;
  vector[N] y_obs;
  vector<lower=0>[N] x_error;
  vector<lower=0>[N] y_error;
}
parameters {
  real alpha;
  real<lower=0> beta;
  real<lower=0> sigma_model;
}
model {
  alpha ~ normal(0, 2);
  beta ~ normal(0, 0.02);
  sigma_model ~ normal(0, 0.8);

  // Integrated Gaussian errors-in-variables likelihood for a linear mean function:
  // if x_true ~ normal(x_obs, x_error), y_obs ~ normal(alpha + beta*x_true, sqrt(y_error^2 + sigma_model^2)),
  // then y_obs | x_obs has variance y_error^2 + sigma_model^2 + beta^2*x_error^2.
  y_obs ~ normal(alpha + beta * x_obs,
                 sqrt(square(y_error) + square(sigma_model) + square(beta * x_error)));
}
generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  vector[N] mu;
  for (n in 1:N) {
    real sd_n = sqrt(square(y_error[n]) + square(sigma_model) + square(beta * x_error[n]));
    mu[n] = alpha + beta * x_obs[n];
    log_lik[n] = normal_lpdf(y_obs[n] | mu[n], sd_n);
    y_rep[n] = normal_rng(mu[n], sd_n);
  }
}
