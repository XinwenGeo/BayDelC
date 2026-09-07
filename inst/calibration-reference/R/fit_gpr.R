# Gaussian Process calibration models --------------------------------------
# Main text GPR: baseline brms GPR, matching the original manuscript logic.
# Supplementary GPR: shrinkage semi-parametric GPR, X-jittered ensemble GPR,
# and custom-Stan latent-input GPR.

standardize_o2_for_gpr <- function(data, center = NULL, scale = NULL) {
  center <- center %||% mean(data$O2, na.rm = TRUE)
  scale <- scale %||% stats::sd(data$O2, na.rm = TRUE)
  data |>
    dplyr::mutate(
      O2_mean_train = center,
      O2_sd_train_scale = scale,
      O2_std = (.data$O2 - center) / scale,
      O2_sd_std = .data$O2_sd / scale
    )
}

make_gpr_formula <- function(model = c("baseline", "shrinkage")) {
  model <- match.arg(model)
  if (model == "baseline") {
    brms::bf(dd13C | se(dd13C_sd, sigma = TRUE) ~ 1 + gp(O2_std, scale = FALSE))
  } else {
    brms::bf(dd13C | se(dd13C_sd, sigma = TRUE) ~ 1 + O2_std + gp(O2_std, scale = FALSE))
  }
}

# Create version-safe priors using brms::get_prior() rather than guessing lscale names.
make_gpr_priors <- function(formula, data, model = c("baseline", "shrinkage"),
                            lscale_prior = "lognormal(log(0.8), 0.4)",
                            sdgp_prior = "exponential(4)",
                            sigma_prior = "exponential(2)",
                            intercept_prior = "normal(0, 2)",
                            b_prior = "normal(0.6, 0.5)",
                            b_positive = FALSE,
                            use_lscale_prior = TRUE) {
  model <- match.arg(model)

  # IMPORTANT for brms >= 2.22.0:
  # Use brms::set_prior() when the prior expression is stored in a character
  # object. brms::prior(sdgp_prior, ...) treats the symbol `sdgp_prior` as the
  # prior itself and leads to errors such as:
  #   sdgp_sdgp_coef[1] ~ sdgp_prior
  #   lscale_lscale_coef[1] ~ lscale_prior
  prior_tbl <- brms::get_prior(
    formula = formula,
    data = data,
    family = stats::gaussian()
  )

  priors <- c(
    brms::set_prior(intercept_prior, class = "Intercept"),
    brms::set_prior(sigma_prior, class = "sigma")
  )

  # GP amplitude. Prefer the coefficient-specific row exposed by get_prior(),
  # e.g. class = "sdgp", coef = "gpO2_std" in brms 2.22.0.
  sdgp_coef <- unique(prior_tbl$coef[
    prior_tbl$class == "sdgp" & !is.na(prior_tbl$coef) & nzchar(prior_tbl$coef)
  ])
  if (length(sdgp_coef) > 0) {
    priors <- c(priors, brms::set_prior(sdgp_prior, class = "sdgp", coef = sdgp_coef[1]))
  } else {
    priors <- c(priors, brms::set_prior(sdgp_prior, class = "sdgp"))
  }

  # GP length-scale. Do not guess names such as lscale_z; use get_prior().
  # In your brms 2.22.0 run, the target is class = "lscale", coef = "gpO2_std".
  if (isTRUE(use_lscale_prior)) {
    lscale_coef <- unique(prior_tbl$coef[
      prior_tbl$class == "lscale" & !is.na(prior_tbl$coef) & nzchar(prior_tbl$coef)
    ])
    if (length(lscale_coef) > 0) {
      priors <- c(priors, brms::set_prior(lscale_prior, class = "lscale", coef = lscale_coef[1]))
    } else {
      has_global_lscale <- any(prior_tbl$class == "lscale" &
                                 (is.na(prior_tbl$coef) | !nzchar(prior_tbl$coef)))
      if (has_global_lscale) {
        priors <- c(priors, brms::set_prior(lscale_prior, class = "lscale"))
      }
    }
  }

  if (model == "shrinkage") {
    # brms 2.22.0 does not allow coef-specific priors together with boundaries
    # (error: "Prior argument 'coef' may not be specified when using boundaries").
    # Default publishable setting: use an informative positive-leaning prior
    # for the O2_std slope without a hard lower bound. If b_positive = TRUE,
    # apply the lower bound at class level, which is valid because O2_std is
    # the only population-level coefficient in the shrinkage formula.
    if (isTRUE(b_positive)) {
      priors <- c(priors, brms::set_prior(b_prior, class = "b", lb = 0))
    } else {
      priors <- c(priors, brms::set_prior(b_prior, class = "b", coef = "O2_std"))
    }
  }

  attr(priors, "prior_table") <- prior_tbl
  priors
}

fit_gpr_brms <- function(data, model_name,
                         gpr_model = c("baseline", "shrinkage"),
                         iter = 4000, warmup = 1000, chains = 4, seed = 1234,
                         adapt_delta = 0.99, max_treedepth = 15,
                         use_lscale_prior = TRUE, b_positive = FALSE,
                         moment_match = TRUE, compute_loo = TRUE) {
  gpr_model <- match.arg(gpr_model)
  ensure_dirs()
  data <- standardize_o2_for_gpr(data)
  formula <- make_gpr_formula(gpr_model)
  priors <- make_gpr_priors(
    formula = formula,
    data = data,
    model = gpr_model,
    use_lscale_prior = use_lscale_prior,
    b_positive = b_positive
  )

  fit <- brms::brm(
    formula = formula,
    data = data,
    family = stats::gaussian(),
    prior = priors,
    iter = iter,
    warmup = warmup,
    chains = chains,
    seed = seed,
    backend = "rstan",
    save_pars = brms::save_pars(all = TRUE),
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  )

  if (isTRUE(compute_loo)) {
    loo_obj <- safe_loo_brms(fit, moment_match = moment_match)
  } else {
    loo_obj <- NULL
  }
  mu <- brms::posterior_epred(fit)
  yhat <- as.numeric(colMeans(mu, na.rm = TRUE))
  yrep <- brms::posterior_predict(fit)

  metrics <- list(
    RMSE = rmse(data$dd13C, yhat),
    Bayesian_R2 = as.numeric(brms::bayes_R2(fit)[1, "Estimate"]),
    PPC_coverage_95 = posterior_coverage(data$dd13C, yrep, 0.95),
    max_pareto_k = extract_pareto_k(loo_obj),
    ELPD = extract_elpd(loo_obj, "Estimate"),
    ELPD_SE = extract_elpd(loo_obj, "SE")
  )

  make_model_bundle(
    model_name = model_name,
    model_type = if (gpr_model == "baseline") "GPR_BASELINE" else "GPR_SHRINKAGE_SEMIPARAMETRIC",
    data = data,
    fit = fit,
    posterior_draws = posterior::as_draws_df(fit),
    log_lik = brms::log_lik(fit),
    y_rep = yrep,
    mu = mu,
    loo = loo_obj,
    metrics = metrics,
    calibration_range = range(data$O2, na.rm = TRUE),
    settings = list(
      brms_version = as.character(utils::packageVersion("brms")),
      formula = deparse(formula$formula),
      gpr_model = gpr_model,
      iter = iter,
      warmup = warmup,
      chains = chains,
      seed = seed,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      compute_loo = compute_loo,
      x_transform = list(mean = unique(data$O2_mean_train), sd = unique(data$O2_sd_train_scale)),
      x_error_model = "not directly included; response-side measurement error included via se()",
      lscale_prior_requested = use_lscale_prior,
      prior_table = attr(priors, "prior_table"),
      prior_summary = brms::prior_summary(fit)
    )
  )
}

fit_gpr_baseline <- function(data, model_name, ...) {
  fit_gpr_brms(data, model_name, gpr_model = "baseline", ...)
}

fit_gpr_shrinkage <- function(data, model_name, ...) {
  fit_gpr_brms(data, model_name, gpr_model = "shrinkage", ...)
}

# Local GPR forward summary used by x-jitter ensembles. This keeps fit_gpr.R
# self-contained even when R/forward.R has not been sourced. If the public
# forward_gpr() exists, use it; otherwise do the same prediction internally.
forward_gpr_safe <- function(bundle, O2_grid, mode = c("predict", "epred"), ndraws = NULL) {
  mode <- match.arg(mode)
  if (exists("forward_gpr", mode = "function", inherits = TRUE)) {
    return(forward_gpr(bundle, O2_grid = O2_grid, mode = mode, ndraws = ndraws))
  }
  if (!inherits(bundle$fit_object, "brmsfit")) stop("forward_gpr_safe requires a brms GPR bundle.")
  tr <- bundle$settings$x_transform
  newdata <- tibble::tibble(
    O2 = as.numeric(O2_grid),
    O2_std = (as.numeric(O2_grid) - tr$mean) / tr$sd,
    dd13C_sd = 0
  )
  pred <- if (mode == "predict") {
    brms::posterior_predict(bundle$fit_object, newdata = newdata, ndraws = ndraws)
  } else {
    brms::posterior_epred(bundle$fit_object, newdata = newdata, ndraws = ndraws)
  }
  tibble::tibble(
    O2 = as.numeric(O2_grid),
    median = apply(pred, 2, stats::median, na.rm = TRUE),
    q025 = apply(pred, 2, stats::quantile, 0.025, na.rm = TRUE),
    q975 = apply(pred, 2, stats::quantile, 0.975, na.rm = TRUE)
  )
}

fit_gpr_xjitter_ensemble <- function(data, model_name,
                                     M = 50,
                                     gpr_model = c("baseline", "shrinkage"),
                                     seed = 3001,
                                     iter = 2000, warmup = 1000, chains = 4,
                                     grid = seq(50, 300, by = 1),
                                     keep_fits = FALSE,
                                     ...) {
  gpr_model <- match.arg(gpr_model)
  ensure_dirs()
  set.seed(seed)
  fits <- vector("list", M)
  curves <- vector("list", M)
  metrics <- vector("list", M)

  for (m in seq_len(M)) {
    dat_m <- data
    dat_m$O2 <- pmax(0, stats::rnorm(nrow(dat_m), mean = dat_m$O2, sd = dat_m$O2_sd))
    # Keep the same response-side errors. O2_sd is represented by the jitter ensemble.
    fit_m <- fit_gpr_brms(
      dat_m,
      model_name = paste0(model_name, "_member", sprintf("%03d", m)),
      gpr_model = gpr_model,
      iter = iter,
      warmup = warmup,
      chains = chains,
      seed = seed + m,
      moment_match = FALSE,
      compute_loo = FALSE,
      ...
    )
    pred_m <- forward_gpr_safe(fit_m, O2_grid = grid, mode = "epred") |>
      dplyr::mutate(member = m)
    curves[[m]] <- pred_m
    metrics[[m]] <- tibble::as_tibble(fit_m$metrics) |> dplyr::mutate(member = m)
    if (keep_fits) fits[[m]] <- fit_m
  }

  curves <- dplyr::bind_rows(curves)
  metrics_df <- dplyr::bind_rows(metrics)
  summary_curve <- curves |>
    dplyr::group_by(.data$O2) |>
    dplyr::summarise(
      median = stats::median(.data$median, na.rm = TRUE),
      q025 = stats::quantile(.data$median, 0.025, na.rm = TRUE),
      q975 = stats::quantile(.data$median, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  bundle <- list(
    model_name = model_name,
    model_type = paste0("GPR_XJITTER_ENSEMBLE_", toupper(gpr_model)),
    bundle_version = baydelc_version(),
    data_version = "O2DD13C_NEW_v0.2",
    training_data = data,
    fits = if (keep_fits) fits else NULL,
    member_curves = curves,
    summary_curve = summary_curve,
    member_metrics = metrics_df,
    settings = list(M = M, gpr_model = gpr_model, seed = seed, iter = iter, warmup = warmup, chains = chains, grid = grid),
    created_time = Sys.time(),
    session_info = sessioninfo::session_info()
  )
  class(bundle) <- c("baydelc_xjitter_bundle", "list")
  bundle
}

fit_gpr_latent_input <- function(data, model_name,
                                 stan_file = "inst/stan/gpr_latent_input.stan",
                                 iter = 4000, warmup = 1000, chains = 4, seed = 5001,
                                 adapt_delta = 0.99, max_treedepth = 15,
                                 compute_loo = TRUE) {
  ensure_dirs()
  dat <- standardize_o2_for_gpr(data)
  stan_data <- list(
    N = nrow(dat),
    x_obs = as.numeric(dat$O2_std),
    x_error = as.numeric(dat$O2_sd_std),
    y_obs = as.numeric(dat$dd13C),
    y_error = as.numeric(dat$dd13C_sd),
    y_bar = mean(dat$dd13C, na.rm = TRUE)
  )
  fit <- rstan::stan(
    file = stan_file,
    data = stan_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    seed = seed,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  )
  draws <- posterior::as_draws_df(fit)
  ext <- rstan::extract(fit)
  loo_obj <- if (isTRUE(compute_loo)) safe_loo_from_loglik(ext$log_lik) else NULL
  yhat <- colMeans(ext$mu, na.rm = TRUE)
  metrics <- list(
    RMSE = rmse(dat$dd13C, yhat),
    Bayesian_R2 = bayes_r2_simple(yhat, dat$dd13C),
    PPC_coverage_95 = posterior_coverage(dat$dd13C, ext$y_rep, 0.95),
    max_pareto_k = if (!is.null(loo_obj)) extract_pareto_k(loo_obj) else NA_real_,
    ELPD = if (!is.null(loo_obj)) extract_elpd(loo_obj, "Estimate") else NA_real_,
    ELPD_SE = if (!is.null(loo_obj)) extract_elpd(loo_obj, "SE") else NA_real_
  )
  make_model_bundle(
    model_name = model_name,
    model_type = "GPR_LATENT_INPUT_STAN",
    data = dat,
    fit = fit,
    posterior_draws = draws,
    log_lik = ext$log_lik,
    y_rep = ext$y_rep,
    mu = ext$mu,
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
      compute_loo = compute_loo,
      x_transform = list(mean = unique(dat$O2_mean_train), sd = unique(dat$O2_sd_train_scale)),
      x_error_model = "latent x_true with x_obs ~ normal(x_true, x_error)"
    ),
    notes = "Computationally intensive sensitivity model; use Sub first."
  )
}
