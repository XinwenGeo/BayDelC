# Bayesian linear calibration models ---------------------------------------
# Main model: Gaussian integrated errors-in-variables BLR.
# Supplementary model: Student-t integrated errors-in-variables BLR.

fit_blr_eiv <- function(data, model_name,
                        stan_file = "inst/stan/blr_eiv_gaussian.stan",
                        model_type = "BLR_EIV_GAUSSIAN",
                        iter = 4000, warmup = 1000, chains = 4, seed = 1234,
                        adapt_delta = 0.98, max_treedepth = 12,
                        save_warmup = FALSE) {
  stopifnot(all(c("O2", "O2_sd", "dd13C", "dd13C_sd") %in% names(data)))
  ensure_dirs()

  stan_data <- list(
    N = nrow(data),
    x_obs = as.numeric(data$O2),
    y_obs = as.numeric(data$dd13C),
    x_error = as.numeric(data$O2_sd),
    y_error = as.numeric(data$dd13C_sd)
  )

  fit <- rstan::stan(
    file = stan_file,
    data = stan_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    seed = seed,
    save_warmup = save_warmup,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  )

  draws <- posterior::as_draws_df(fit)
  ext <- rstan::extract(fit)
  log_lik <- ext$log_lik
  y_rep <- ext$y_rep
  mu <- ext$mu
  loo_obj <- safe_loo_from_loglik(log_lik)

  yhat <- colMeans(mu, na.rm = TRUE)
  metrics <- list(
    RMSE = rmse(data$dd13C, yhat),
    Bayesian_R2 = bayes_r2_simple(yhat, data$dd13C),
    PPC_coverage_95 = posterior_coverage(data$dd13C, y_rep, 0.95),
    max_pareto_k = extract_pareto_k(loo_obj),
    ELPD = extract_elpd(loo_obj, "Estimate"),
    ELPD_SE = extract_elpd(loo_obj, "SE")
  )

  make_model_bundle(
    model_name = model_name,
    model_type = model_type,
    data = data,
    fit = fit,
    posterior_draws = draws,
    log_lik = log_lik,
    y_rep = y_rep,
    mu = mu,
    loo = loo_obj,
    metrics = metrics,
    calibration_range = range(data$O2, na.rm = TRUE),
    settings = list(
      stan_file = stan_file,
      iter = iter,
      warmup = warmup,
      chains = chains,
      seed = seed,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      likelihood = model_type,
      x_error_model = "integrated linear errors-in-variables: y variance includes beta^2 * x_error^2"
    )
  )
}

fit_blr_eiv_gaussian <- function(data, model_name, ...) {
  fit_blr_eiv(
    data = data,
    model_name = model_name,
    stan_file = "inst/stan/blr_eiv_gaussian.stan",
    model_type = "BLR_EIV_GAUSSIAN",
    ...
  )
}

fit_blr_eiv_student_t <- function(data, model_name, ...) {
  fit_blr_eiv(
    data = data,
    model_name = model_name,
    stan_file = "inst/stan/blr_eiv_student_t.stan",
    model_type = "BLR_EIV_STUDENT_T",
    adapt_delta = 0.99,
    max_treedepth = 15,
    ...
  )
}
