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
  real<lower=2> nu;
}
model {
  alpha ~ normal(0, 2);
  beta ~ normal(0, 0.02);
  sigma_model ~ normal(0, 0.8);
  nu ~ gamma(2, 0.1);

  y_obs ~ student_t(nu, alpha + beta * x_obs,
                    sqrt(square(y_error) + square(sigma_model) + square(beta * x_error)));
}
generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  vector[N] mu;
  for (n in 1:N) {
    real sd_n = sqrt(square(y_error[n]) + square(sigma_model) + square(beta * x_error[n]));
    mu[n] = alpha + beta * x_obs[n];
    log_lik[n] = student_t_lpdf(y_obs[n] | nu, mu[n], sd_n);
    y_rep[n] = student_t_rng(nu, mu[n], sd_n);
  }
}
