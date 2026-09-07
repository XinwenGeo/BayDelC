# BayDelC v0.3 unified inverse interface -----------------------------------
# The helper functions in forward_psm.R are loaded first through DESCRIPTION's
# Collate field.  Keeping package files free of source() calls makes BayDelC
# work consistently from source, binary, and installed packages.

#' Describe uncertainty components for a Delta-delta-13-C observation
#'
#' This helper combines epifaunal and infaunal analytical or sampling
#' uncertainties, their correlation, and optional pairing or other terms. All
#' values are interpreted as one-standard-deviation errors unless `scale =
#' "se"`; the current implementation treats supplied standard errors as the
#' standard deviations of the corresponding observation-level estimates.
#'
#' @param epi_total_sd,infa_total_sd Total standard deviations for the
#'   epifaunal and infaunal measurements. Do not combine these with their
#'   corresponding component arguments.
#' @param epi_measurement_sd,epi_sample_sd Epifaunal measurement and sampling
#'   standard deviations.
#' @param infa_measurement_sd,infa_sample_sd Infaunal measurement and sampling
#'   standard deviations.
#' @param pairing_sd,other_sd Additional independent standard-deviation terms.
#' @param rho_epi_infa Correlation between epifaunal and infaunal errors.
#' @param scale Input scale, currently `"sd"` or `"se"`.
#' @return An object suitable for `inverse_psm(error_components = ...)`.
#' @export
psm_error_components <- function(
    epi_total_sd = NULL,
    infa_total_sd = NULL,
    epi_measurement_sd = NULL,
    epi_sample_sd = NULL,
    infa_measurement_sd = NULL,
    infa_sample_sd = NULL,
    pairing_sd = 0,
    other_sd = 0,
    rho_epi_infa = 0,
    scale = c("sd", "se")) {
  scale <- match.arg(scale)
  if (!is.null(epi_total_sd) && (!is.null(epi_measurement_sd) || !is.null(epi_sample_sd))) {
    stop("Use either `epi_total_sd` or epifaunal measurement/sample components, not both.", call. = FALSE)
  }
  if (!is.null(infa_total_sd) && (!is.null(infa_measurement_sd) || !is.null(infa_sample_sd))) {
    stop("Use either `infa_total_sd` or infaunal measurement/sample components, not both.", call. = FALSE)
  }
  structure(list(
    epi_total_sd = epi_total_sd,
    infa_total_sd = infa_total_sd,
    epi_measurement_sd = epi_measurement_sd,
    epi_sample_sd = epi_sample_sd,
    infa_measurement_sd = infa_measurement_sd,
    infa_sample_sd = infa_sample_sd,
    pairing_sd = pairing_sd,
    other_sd = other_sd,
    rho_epi_infa = rho_epi_infa,
    scale = scale
  ), class = c("psm_error_components", "list"))
}

.psm_component_vector <- function(x, n, name, default = 0) {
  if (is.null(x)) return(rep(default, n))
  x <- .psm_recycle_numeric(x, n, name)
  if (any(!is.finite(x) | x < 0)) stop(sprintf("`%s` must be finite and non-negative.", name), call. = FALSE)
  x
}

.psm_components_to_sd <- function(components, n) {
  if (is.data.frame(components)) components <- as.list(components)
  if (!is.list(components)) stop("`error_components` must be a list/data frame or output of `psm_error_components()`.", call. = FALSE)
  if (!is.null(components$epi_total_sd) && (!is.null(components$epi_measurement_sd) || !is.null(components$epi_sample_sd))) stop("Epifaunal total SD conflicts with component SDs.", call. = FALSE)
  if (!is.null(components$infa_total_sd) && (!is.null(components$infa_measurement_sd) || !is.null(components$infa_sample_sd))) stop("Infaunal total SD conflicts with component SDs.", call. = FALSE)
  epi <- if (!is.null(components$epi_total_sd)) {
    .psm_component_vector(components$epi_total_sd, n, "epi_total_sd")
  } else {
    sqrt(.psm_component_vector(components$epi_measurement_sd, n, "epi_measurement_sd")^2 +
           .psm_component_vector(components$epi_sample_sd, n, "epi_sample_sd")^2)
  }
  infa <- if (!is.null(components$infa_total_sd)) {
    .psm_component_vector(components$infa_total_sd, n, "infa_total_sd")
  } else {
    sqrt(.psm_component_vector(components$infa_measurement_sd, n, "infa_measurement_sd")^2 +
           .psm_component_vector(components$infa_sample_sd, n, "infa_sample_sd")^2)
  }
  rho <- .psm_recycle_numeric(components$rho_epi_infa %||% 0, n, "rho_epi_infa")
  if (any(!is.finite(rho) | rho < -1 | rho > 1)) stop("`rho_epi_infa` must lie between -1 and 1.", call. = FALSE)
  pairing <- .psm_component_vector(components$pairing_sd, n, "pairing_sd")
  other <- .psm_component_vector(components$other_sd, n, "other_sd")
  variance <- epi^2 + infa^2 - 2 * rho * epi * infa + pairing^2 + other^2
  if (any(variance < -1e-12)) stop("The component covariance specification produces a negative variance.", call. = FALSE)
  sqrt(pmax(0, variance))
}

.psm_resolve_observation_error <- function(
    n_obs, dd13c_sd, glob_site_sd, error_components, dd13c_draws,
    error_method, validate, draw_weights = NULL) {
  draw_mat <- .psm_as_draw_matrix(dd13c_draws, n_obs, "dd13c_draws", validate)
  if (!is.null(dd13c_sd)) dd13c_sd <- .psm_recycle_numeric(dd13c_sd, n_obs, "dd13c_sd")
  if (!is.null(glob_site_sd)) glob_site_sd <- .psm_recycle_numeric(glob_site_sd, n_obs, "glob_site_sd")
  reps <- c(provided_total = !is.null(dd13c_sd), lu_site = !is.null(glob_site_sd), components = !is.null(error_components), draws = !is.null(draw_mat))
  if (error_method == "auto") {
    if (sum(reps) > 1L) .psm_issue("Conflicting Delta-delta13C uncertainty representations were supplied.", validate)
    error_method <- if (!is.null(draw_mat)) "draws" else if (!is.null(dd13c_sd)) "provided_total" else if (!is.null(error_components)) "components" else if (!is.null(glob_site_sd)) "lu_site" else "h15_proxy_default"
  }
  required <- switch(error_method,
    provided_total = !is.null(dd13c_sd),
    lu_site = !is.null(glob_site_sd),
    components = !is.null(error_components),
    draws = !is.null(draw_mat),
    h15_proxy_default = TRUE,
    FALSE
  )
  if (!required) stop(sprintf("Required uncertainty input is missing for `error_method = '%s'`.", error_method), call. = FALSE)
  irrelevant <- switch(error_method,
    provided_total = reps[c("lu_site", "components", "draws")],
    lu_site = reps[c("provided_total", "components", "draws")],
    components = reps[c("provided_total", "lu_site", "draws")],
    draws = reps[c("provided_total", "lu_site", "components")],
    h15_proxy_default = reps,
    logical()
  )
  if (any(irrelevant)) .psm_issue("Inputs incompatible with the selected `error_method` were also supplied.", validate)
  sy <- switch(error_method,
    provided_total = dd13c_sd,
    lu_site = sqrt(0.05^2 + glob_site_sd^2),
    components = .psm_components_to_sd(error_components, n_obs),
    h15_proxy_default = rep(sqrt(0.05^2 + 0.05^2 + 0.05^2), n_obs),
    draws = rep(NA_real_, n_obs)
  )
  if (!is.null(sy) && any(!is.na(sy) & (!is.finite(sy) | sy < 0))) stop("Resolved Delta-delta13C SD values must be finite and non-negative.", call. = FALSE)
  w <- if (!is.null(draw_mat)) .psm_normalize_weights(draw_weights, ncol(draw_mat), "dd13c_draw_weights") else NULL
  list(method = error_method, sd = sy, draws = draw_mat, weights = w)
}

.psm_logmeanexp_rows <- function(x) {
  m <- apply(x, 2L, max)
  m + log(colMeans(exp(sweep(x, 2L, m, "-"))))
}

.psm_log_density <- function(y, mu, scale, student = FALSE, nu = NULL) {
  scale <- pmax(as.numeric(scale), 1e-9)
  if (student) {
    nu <- ifelse(is.finite(nu) & nu > 2, nu, 5)
    stats::dt((y - mu) / scale, df = nu, log = TRUE) - log(scale)
  } else {
    stats::dnorm(y, mean = mu, sd = scale, log = TRUE)
  }
}

.psm_grid_cell_widths <- function(x) {
  if (length(x) == 1L) return(1)
  d <- diff(x)
  c(d[1L] / 2, (d[-1L] + d[-length(d)]) / 2, d[length(d)] / 2)
}

.psm_density_to_mass <- function(density, grid) {
  density[!is.finite(density) | density < 0] <- 0
  mass <- density * .psm_grid_cell_widths(grid)
  if (sum(mass) <= 0) stop("The posterior has zero numerical mass on the O2 grid.", call. = FALSE)
  mass / sum(mass)
}

.psm_quantile_mass <- function(grid, mass, probs = c(0.025, 0.5, 0.975)) {
  cdf <- cumsum(mass)
  vapply(probs, function(p) stats::approx(c(0, cdf), c(grid[1L], grid), xout = p, ties = "ordered", rule = 2)$y, numeric(1))
}

.psm_weighted_quantile <- function(x, w, probs = c(0.025, 0.5, 0.975)) {
  ok <- is.finite(x) & is.finite(w) & w >= 0
  x <- x[ok]; w <- w[ok]
  if (!length(x) || sum(w) <= 0) return(rep(NA_real_, length(probs)))
  o <- order(x); x <- x[o]; w <- w[o] / sum(w[o])
  cdf <- cumsum(w)
  vapply(probs, function(p) stats::approx(c(0, cdf), c(x[1L], x), xout = p, ties = "ordered", rule = 2)$y, numeric(1))
}

.psm_prior_normal_args <- function(prior_args, n_obs) {
  if (is.null(prior_args$mean) || is.null(prior_args$sd)) stop("Normal prior requires `prior_args$mean` and `prior_args$sd`.", call. = FALSE)
  mu <- .psm_recycle_numeric(prior_args$mean, n_obs, "prior_args$mean")
  sd <- .psm_recycle_numeric(prior_args$sd, n_obs, "prior_args$sd")
  if (any(!is.finite(sd) | sd <= 0)) stop("Normal prior SD must be positive.", call. = FALSE)
  bounds <- prior_args$bounds
  if (!is.null(bounds)) bounds <- .psm_validate_bounds(bounds, "prior_args$bounds")
  list(mean = mu, sd = sd, bounds = bounds)
}

.psm_build_grid <- function(prior, prior_args, o2_bounds, o2_step, o2_grid, y, sy, model, algorithm) {
  if (!is.null(o2_grid)) {
    grid <- sort(unique(as.numeric(o2_grid)))
    if (length(grid) < 2L || any(!is.finite(grid))) stop("`o2_grid` must contain at least two finite unique values.", call. = FALSE)
    return(grid)
  }
  if (!is.finite(o2_step) || o2_step <= 0) stop("`o2_step` must be positive.", call. = FALSE)
  bounds <- .psm_validate_bounds(o2_bounds, "o2_bounds")
  if (prior == "normal") {
    pa <- .psm_prior_normal_args(prior_args, length(y))
    if (!is.null(pa$bounds)) {
      bounds <- pa$bounds
    } else {
      bounds <- c(min(pa$mean - 5 * pa$sd), max(pa$mean + 5 * pa$sd))
      if (algorithm == "BLR") {
        pars <- .psm_blr_parameter_draws(model, 1000, 7401)
        a <- stats::median(pars$alpha); b <- stats::median(pars$beta); s <- stats::median(pars$sigma_model)
        if (is.finite(b) && abs(b) > 1e-8) {
          yc <- (y - a) / b
          spread <- 6 * sqrt(ifelse(is.finite(sy), sy^2, 0) + s^2) / abs(b)
          bounds <- range(c(bounds, yc - spread, yc + spread), finite = TRUE)
        }
      }
    }
  } else if (prior == "density" && !is.null(prior_args$x)) {
    bounds <- range(as.numeric(prior_args$x), finite = TRUE)
  }
  seq(floor(bounds[1] / o2_step) * o2_step, ceiling(bounds[2] / o2_step) * o2_step, by = o2_step)
}

.psm_prior_density <- function(prior, prior_args, grid, row, n_obs, data = NULL) {
  if (prior == "none") return(rep(1, length(grid)))
  if (prior == "uniform_physical") return(ifelse(grid >= 0 & grid <= 300, 1, 0))
  if (prior == "normal") {
    pa <- .psm_prior_normal_args(prior_args, n_obs)
    d <- stats::dnorm(grid, pa$mean[row], pa$sd[row])
    if (!is.null(pa$bounds)) d[grid < pa$bounds[1] | grid > pa$bounds[2]] <- 0
    return(d)
  }
  if (prior == "density") {
    obj <- prior_args$density %||% prior_args
    if (is.function(obj)) {
      ans <- try(obj(grid, row = row, data = data), silent = TRUE)
      if (inherits(ans, "try-error")) ans <- obj(grid)
      return(as.numeric(ans))
    }
    if (is.data.frame(obj)) {
      xcol <- intersect(c("O2", "o2", "x"), names(obj)); dcol <- intersect(c("density", "probability", "weight", "y"), names(obj))
      if (!length(xcol) || !length(dcol)) stop("Prior density data frame needs O2/x and density columns.", call. = FALSE)
      return(stats::approx(obj[[xcol[1L]]], obj[[dcol[1L]]], xout = grid, yleft = 0, yright = 0)$y)
    }
    obj <- as.numeric(obj)
    if (length(obj) != length(grid)) stop("Numeric prior density must match the O2 grid length.", call. = FALSE)
    return(obj)
  }
  if (prior == "custom") {
    fun <- prior_args$fun %||% prior_args$density
    if (!is.function(fun)) stop("Custom prior requires `prior_args$fun`.", call. = FALSE)
    ans <- try(fun(grid = grid, row = row, data = data), silent = TRUE)
    if (inherits(ans, "try-error")) ans <- fun(grid)
    return(as.numeric(ans))
  }
  stop(sprintf("Unsupported prior `%s`.", prior), call. = FALSE)
}

.psm_gpr_mean_at_points <- function(model, x, draw_ids) {
  .psm_gpr_mean_grid(model, x, draw_ids)
}

.psm_inverse_ensemble_prior <- function(
    y, obs, prior_args, model, algorithm, sites, include_model_residual,
    n_draw, integration_draws, seed, return_draws, return_density, calib_range) {
  prior_mat <- .psm_as_draw_matrix(prior_args$draws, length(y), "prior_args$draws", "strict")
  prior_w <- .psm_normalize_weights(prior_args$weights %||% NULL, ncol(prior_mat), "prior_args$weights")
  s <- min(as.integer(integration_draws), max(1L, n_draw))
  if (algorithm == "BLR") {
    pars <- .psm_blr_parameter_draws(model, s, seed + 11L)
    student <- .psm_is_student_t(model, pars)
    draw_ids <- NULL
  } else {
    n_avail <- .psm_gpr_available_draws(model)
    draw_ids <- .psm_sample_rows(n_avail, s, seed + 11L)
    sigma <- if (include_model_residual) .psm_gpr_sigma_draws(model, draw_ids) else rep(0, s)
    student <- FALSE
  }
  obs_idx <- if (!is.null(obs$draws)) sample.int(ncol(obs$draws), s, replace = s > ncol(obs$draws), prob = obs$weights) else NULL
  log_like <- matrix(0, nrow = s, ncol = ncol(prior_mat))
  for (i in seq_along(y)) {
    yi <- if (is.null(obs_idx)) rep(y[i], s) else obs$draws[i, obs_idx]
    syi <- if (is.null(obs_idx)) obs$sd[i] else 0
    if (algorithm == "BLR") {
      mu <- outer(pars$beta, prior_mat[i, ], "*") + pars$alpha
      scale <- sqrt((if (include_model_residual) pars$sigma_model else 0)^2 + syi^2)
      ll <- matrix(NA_real_, nrow = s, ncol = ncol(prior_mat))
      for (k in seq_len(s)) ll[k, ] <- .psm_log_density(yi[k], mu[k, ], scale[k], student, pars$nu[k])
    } else {
      mu <- .psm_gpr_mean_at_points(model, prior_mat[i, ], draw_ids)
      scale <- sqrt(sigma^2 + syi^2)
      ll <- matrix(NA_real_, nrow = s, ncol = ncol(prior_mat))
      for (k in seq_len(s)) ll[k, ] <- .psm_log_density(yi[k], mu[k, ], scale[k], FALSE, NA)
    }
    log_like <- log_like + ll
  }
  log_w <- log(prior_w + 1e-300) + .psm_logmeanexp_rows(log_like)
  log_w <- log_w - max(log_w)
  post_w <- exp(log_w); post_w <- post_w / sum(post_w)
  summary <- do.call(rbind, lapply(seq_along(y), function(i) {
    q <- .psm_weighted_quantile(prior_mat[i, ], post_w)
    data.frame(row_id = i, Site = sites[i], input_value = y[i], input_sd = obs$sd[i],
               mean = sum(prior_mat[i, ] * post_w), median = q[2L], q025 = q[1L], q975 = q[3L], stringsAsFactors = FALSE)
  }))
  set.seed(seed + 12L)
  idx <- sample.int(ncol(prior_mat), n_draw, replace = TRUE, prob = post_w)
  out_draws <- prior_mat[, idx, drop = FALSE]
  density <- if (return_density) do.call(rbind, lapply(seq_along(y), function(i) data.frame(row_id = i, Site = sites[i], ensemble_id = seq_len(ncol(prior_mat)), O2 = prior_mat[i, ], density = post_w))) else NULL
  weighted_indicator <- function(condition) rowSums(condition * matrix(post_w, nrow = nrow(prior_mat), ncol = ncol(prior_mat), byrow = TRUE))
  probabilities <- data.frame(
    P_below_zero = weighted_indicator(prior_mat < 0),
    P_above_300 = weighted_indicator(prior_mat > 300),
    P_below_50 = weighted_indicator(prior_mat < 50),
    P_outside_calibration = weighted_indicator(prior_mat < calib_range[1] | prior_mat > calib_range[2])
  )
  list(summary = summary, draws = if (return_draws) out_draws else NULL, density = density,
       probabilities = probabilities, posterior_weights = post_w, effective_sample_size = 1 / sum(post_w^2))
}

#' Reconstruct bottom-water oxygen from Delta-delta-13-C
#'
#' `inverse_psm()` propagates calibration and observation uncertainty to
#' estimate bottom-water oxygen (BWO). BLR can use direct analytic inversion;
#' BLR and GPR can both use a grid-based Bayesian inversion with an optional
#' prior.
#'
#' @param data Optional data frame containing Delta-delta-13-C observations and
#'   optional uncertainty and site columns.
#' @param dd13c Numeric Delta-delta-13-C observations in per mille.
#' @param dd13c_sd Total one-standard-deviation uncertainty for `dd13c`.
#' @param glob_site_sd Site-specific infaunal standard deviations used by the
#'   Lu-style error representation.
#' @param error_components Output of `psm_error_components()` or a compatible
#'   list/data frame.
#' @param dd13c_draws Optional matrix of observation draws, with observations
#'   in rows and realizations in columns.
#' @param dd13c_draw_weights Optional weights for columns of `dd13c_draws`.
#' @param dd13c_col,dd13c_sd_col,glob_site_sd_col,site_col Column names used
#'   when `data` is supplied.
#' @param error_method Observation-error representation. `"auto"` selects from
#'   the supplied inputs.
#' @param algorithm Calibration family: `"BLR"` or `"GPR"`.
#' @param calibration Calibration dataset: restricted-taxon `"Sub"`,
#'   mixed-taxon `"All"`, or `"custom"`.
#' @param model Optional custom model list or path to an RDS model bundle.
#' @param inversion `"strict"` for analytic BLR inversion, `"grid"` for
#'   Bayesian grid inversion, or `"auto"`.
#' @param prior Prior on BWO. See the package vignette for the accepted
#'   `prior_args` structures.
#' @param prior_args Named list configuring the selected prior.
#' @param o2_bounds,o2_step,o2_grid Bounds, spacing, or an explicit numerical
#'   grid for grid-based inversion, in micromoles per kilogram.
#' @param dd13c_convention Sign convention. This version supports only
#'   `"epifaunal_minus_infaunal"`.
#' @param include_model_residual Include residual calibration variability.
#' @param n_draw Number of output Monte Carlo draws.
#' @param integration_draws Number of posterior model draws used for numerical
#'   likelihood integration.
#' @param seed Random seed.
#' @param return_draws,return_density Return posterior draws or grid density.
#' @param validate Validation policy: `"strict"`, `"warn"`, or `"none"`.
#' @return An object of class `inverse_psm_result` containing `summary`, optional
#'   draws/densities, diagnostics, and model metadata.
#' @export
#' @examples
#' inverse_psm(dd13c = c(1.10, 1.65), dd13c_sd = 0.08,
#'             algorithm = "BLR", n_draw = 500)
#' inverse_psm(dd13c = 1.50, dd13c_sd = 0.08, algorithm = "GPR",
#'             inversion = "grid", prior = "uniform_physical",
#'             integration_draws = 100)
inverse_psm <- function(
    data = NULL,
    dd13c = NULL,
    dd13c_sd = NULL,
    glob_site_sd = NULL,
    error_components = NULL,
    dd13c_draws = NULL,
    dd13c_draw_weights = NULL,
    dd13c_col = "dd13C",
    dd13c_sd_col = NULL,
    glob_site_sd_col = NULL,
    site_col = "Site",
    error_method = c("auto", "provided_total", "h15_proxy_default", "lu_site", "components", "draws"),
    algorithm = c("BLR", "GPR"),
    calibration = c("Sub", "All", "custom"),
    model = NULL,
    inversion = c("auto", "strict", "grid"),
    prior = c("none", "uniform_physical", "normal", "draws", "density", "custom"),
    prior_args = NULL,
    o2_bounds = c(0, 300),
    o2_step = 1,
    o2_grid = NULL,
    dd13c_convention = "epifaunal_minus_infaunal",
    include_model_residual = TRUE,
    n_draw = 2000,
    integration_draws = NULL,
    seed = 123,
    return_draws = FALSE,
    return_density = FALSE,
    validate = c("strict", "warn", "none")) {

  error_method <- match.arg(error_method)
  algorithm <- match.arg(algorithm)
  calibration <- match.arg(calibration)
  inversion <- match.arg(inversion)
  prior <- match.arg(prior)
  validate <- match.arg(validate)
  if (!identical(dd13c_convention, "epifaunal_minus_infaunal")) stop("Only `epifaunal_minus_infaunal` is supported in this version.", call. = FALSE)
  n_draw <- as.integer(n_draw)
  if (!is.finite(n_draw) || n_draw < 1L) stop("`n_draw` must be a positive integer.", call. = FALSE)
  integration_draws <- as.integer(integration_draws %||% min(500L, n_draw))
  if (!is.finite(integration_draws) || integration_draws < 1L) stop("`integration_draws` must be positive.", call. = FALSE)
  prior_args <- prior_args %||% list()

  if (!is.null(data)) {
    data <- as.data.frame(data)
    if (!is.null(dd13c)) .psm_issue("Supply Delta-delta13C through either `data` or `dd13c`, not both.", validate)
    yc <- .psm_detect_column(data, dd13c_col, c("dd13C", "dd13c", "Ddl3C"), TRUE, "Delta-delta13C column")
    dd13c <- as.numeric(data[[yc]])
    if (!is.null(dd13c_sd_col)) {
      sc <- .psm_detect_column(data, dd13c_sd_col, c("dd13C_sd", "dd13c_sd"), TRUE, "Delta-delta13C SD column")
      if (!is.null(dd13c_sd)) .psm_issue("Supply Delta-delta13C SD through either `data` or `dd13c_sd`, not both.", validate)
      dd13c_sd <- as.numeric(data[[sc]])
    } else if (error_method == "auto" && is.null(dd13c_sd)) {
      auto_sc <- .psm_detect_column(data, NULL, c("dd13C_sd", "dd13c_sd"), FALSE, "Delta-delta13C SD column")
      if (!is.null(auto_sc)) dd13c_sd <- as.numeric(data[[auto_sc]])
    }
    if (!is.null(glob_site_sd_col)) {
      gc <- .psm_detect_column(data, glob_site_sd_col, c("glob_site_sd", "Glob_site_sd"), TRUE, "Globobulimina site SD column")
      if (!is.null(glob_site_sd)) .psm_issue("Supply Globobulimina site SD through either `data` or `glob_site_sd`, not both.", validate)
      glob_site_sd <- as.numeric(data[[gc]])
    }
  }
  dd13c <- as.numeric(dd13c)
  if (!length(dd13c) || any(!is.finite(dd13c))) stop("`dd13c` must contain finite values.", call. = FALSE)
  n_obs <- length(dd13c)
  sites <- .psm_site_labels(data, site_col, n_obs)
  obs <- .psm_resolve_observation_error(n_obs, dd13c_sd, glob_site_sd, error_components, dd13c_draws, error_method, validate, dd13c_draw_weights)

  resolved <- .psm_resolve_model(model, algorithm, calibration)
  bundle <- resolved$bundle
  algorithm_used <- resolved$algorithm
  calib_range <- .psm_calibration_range(bundle)
  if (inversion == "auto") inversion <- if (algorithm_used == "BLR" && prior == "none") "strict" else "grid"
  if (inversion == "strict" && algorithm_used != "BLR") stop("Strict analytic inversion is available only for BLR.", call. = FALSE)
  if (inversion == "strict" && prior != "none") stop("Strict analytic inversion is compatible only with `prior = 'none'`; use grid inversion for a prior.", call. = FALSE)
  if (prior == "draws" && is.null(prior_args$draws)) stop("`prior = 'draws'` requires `prior_args$draws`.", call. = FALSE)

  if (inversion == "strict") {
    pars <- .psm_blr_parameter_draws(bundle, n_draw, seed + 1L)
    set.seed(seed + 2L)
    if (!is.null(obs$draws)) {
      idx <- sample.int(ncol(obs$draws), n_draw, replace = n_draw > ncol(obs$draws), prob = obs$weights)
      y_draw <- obs$draws[, idx, drop = FALSE]
    } else {
      y_draw <- matrix(stats::rnorm(n_obs * n_draw, mean = rep(dd13c, n_draw), sd = rep(obs$sd, n_draw)), nrow = n_obs)
    }
    eps <- matrix(0, nrow = n_obs, ncol = n_draw)
    if (include_model_residual) {
      if (.psm_is_student_t(bundle, pars)) {
        nu <- ifelse(is.finite(pars$nu) & pars$nu > 2, pars$nu, 5)
        eps <- matrix(stats::rt(n_obs * n_draw, df = rep(nu, each = n_obs)), nrow = n_obs)
        eps <- sweep(eps, 2L, pars$sigma_model, "*")
      } else {
        eps <- matrix(stats::rnorm(n_obs * n_draw), nrow = n_obs)
        eps <- sweep(eps, 2L, pars$sigma_model, "*")
      }
    }
    denom <- pars$beta
    invalid_beta <- !is.finite(denom) | abs(denom) < 1e-10
    denom[invalid_beta] <- NA_real_
    o2_draws <- sweep(y_draw - eps, 2L, pars$alpha, "-")
    o2_draws <- sweep(o2_draws, 2L, denom, "/")
    summary <- .psm_summary_matrix(o2_draws, sites, dd13c, obs$sd)
    prior_requested <- "none"; prior_effective <- "none"
    density <- NULL
    diag_below0 <- rowMeans(o2_draws < 0, na.rm = TRUE)
    diag_above300 <- rowMeans(o2_draws > 300, na.rm = TRUE)
    diag_below50 <- rowMeans(o2_draws < 50, na.rm = TRUE)
    diag_outside <- rowMeans(o2_draws < calib_range[1] | o2_draws > calib_range[2], na.rm = TRUE)
    extra_diag <- list(P_beta_nonpositive = mean(pars$beta <= 0, na.rm = TRUE), P_beta_invalid = mean(invalid_beta))
  } else if (prior == "draws") {
    ens <- .psm_inverse_ensemble_prior(dd13c, obs, prior_args, bundle, algorithm_used, sites,
                                       include_model_residual, n_draw, integration_draws, seed,
                                       return_draws, return_density, calib_range)
    summary <- ens$summary
    o2_draws <- ens$draws
    density <- ens$density
    prior_requested <- "draws"; prior_effective <- "ensemble_prior"
    diag_below0 <- ens$probabilities$P_below_zero
    diag_above300 <- ens$probabilities$P_above_300
    diag_below50 <- ens$probabilities$P_below_50
    diag_outside <- ens$probabilities$P_outside_calibration
    extra_diag <- list(prior_effective_sample_size = ens$effective_sample_size)
  } else {
    grid <- .psm_build_grid(prior, prior_args, o2_bounds, o2_step, o2_grid, dd13c, obs$sd, bundle, algorithm_used)
    s <- integration_draws
    set.seed(seed + 3L)
    obs_idx <- if (!is.null(obs$draws)) sample.int(ncol(obs$draws), s, replace = s > ncol(obs$draws), prob = obs$weights) else NULL
    if (algorithm_used == "BLR") {
      pars <- .psm_blr_parameter_draws(bundle, s, seed + 4L)
      student <- .psm_is_student_t(bundle, pars)
      mu_grid <- outer(pars$beta, grid, "*") + pars$alpha
      sigma_model <- if (include_model_residual) pars$sigma_model else rep(0, s)
    } else {
      n_avail <- .psm_gpr_available_draws(bundle)
      draw_ids <- .psm_sample_rows(n_avail, s, seed + 4L)
      mu_grid <- .psm_gpr_mean_grid(bundle, grid, draw_ids)
      sigma_model <- if (include_model_residual) .psm_gpr_sigma_draws(bundle, draw_ids) else rep(0, s)
      student <- FALSE; pars <- NULL
    }
    summaries <- vector("list", n_obs); densities <- if (return_density) vector("list", n_obs) else NULL
    o2_draws <- if (return_draws) matrix(NA_real_, nrow = n_obs, ncol = n_draw) else NULL
    diag_below0 <- diag_above300 <- diag_below50 <- diag_outside <- numeric(n_obs)
    for (i in seq_len(n_obs)) {
      yi <- if (is.null(obs_idx)) rep(dd13c[i], s) else obs$draws[i, obs_idx]
      syi <- if (is.null(obs_idx)) obs$sd[i] else 0
      loglik <- matrix(NA_real_, nrow = s, ncol = length(grid))
      scale <- sqrt(sigma_model^2 + syi^2)
      for (k in seq_len(s)) {
        loglik[k, ] <- .psm_log_density(yi[k], mu_grid[k, ], scale[k], student, if (is.null(pars)) NA else pars$nu[k])
      }
      log_like <- .psm_logmeanexp_rows(loglik)
      like <- exp(log_like - max(log_like, na.rm = TRUE))
      pr <- .psm_prior_density(prior, prior_args, grid, i, n_obs, data)
      mass <- .psm_density_to_mass(like * pr, grid)
      q <- .psm_quantile_mass(grid, mass)
      summaries[[i]] <- data.frame(row_id = i, Site = sites[i], input_value = dd13c[i], input_sd = obs$sd[i],
                                    mean = sum(grid * mass), median = q[2L], q025 = q[1L], q975 = q[3L], stringsAsFactors = FALSE)
      if (return_density) densities[[i]] <- data.frame(row_id = i, Site = sites[i], O2 = grid, density = mass / .psm_grid_cell_widths(grid))
      diag_below0[i] <- sum(mass[grid < 0])
      diag_above300[i] <- sum(mass[grid > 300])
      diag_below50[i] <- sum(mass[grid < 50])
      diag_outside[i] <- sum(mass[grid < calib_range[1] | grid > calib_range[2]])
      if (return_draws) o2_draws[i, ] <- sample(grid, n_draw, replace = TRUE, prob = mass)
    }
    summary <- do.call(rbind, summaries)
    density <- if (return_density) do.call(rbind, densities) else NULL
    prior_requested <- prior
    prior_effective <- if (prior == "none") "bounded_flat_reference" else if (prior == "normal" && is.null(prior_args$bounds)) "unbounded_normal_on_adaptive_grid" else if (prior == "normal") "truncated_normal" else prior
    extra_diag <- list(grid_min = min(grid), grid_max = max(grid), o2_step = stats::median(diff(grid)))
  }

  if (!exists("diag_below0")) {
    diag_below0 <- diag_above300 <- diag_below50 <- diag_outside <- rep(NA_real_, n_obs)
  }
  summary$dd13c <- dd13c
  summary$dd13c_sd <- obs$sd
  summary$algorithm <- algorithm_used
  summary$calibration <- calibration
  summary$inversion <- inversion
  summary$prior_requested <- prior_requested
  summary$prior_effective <- prior_effective
  summary$error_method <- obs$method
  summary$P_below_zero <- diag_below0
  summary$P_above_300 <- diag_above300
  summary$P_below_50 <- diag_below50
  summary$P_outside_calibration <- diag_outside

  diagnostics <- summary[, c("row_id", "Site", "P_below_zero", "P_above_300", "P_below_50", "P_outside_calibration"), drop = FALSE]
  input_out <- if (is.null(data)) data.frame(row_id = seq_len(n_obs), Site = sites, dd13c = dd13c, stringsAsFactors = FALSE) else data
  result <- list(
    input = input_out,
    summary = summary,
    draws = if (return_draws) o2_draws else NULL,
    density = density,
    metadata = c(list(
      direction = "inverse",
      model_name = bundle$model_name %||% "custom",
      model_type = bundle$model_type %||% algorithm_used,
      algorithm = algorithm_used,
      calibration = calibration,
      calibration_range = calib_range,
      inversion = inversion,
      prior_requested = prior_requested,
      prior_effective = prior_effective,
      error_method = obs$method,
      include_model_residual = include_model_residual,
      dd13c_convention = dd13c_convention,
      draw_orientation = "rows = observations; columns = joint realizations",
      n_draw = n_draw,
      integration_draws = integration_draws,
      seed = seed
    ), extra_diag),
    diagnostics = diagnostics
  )
  class(result) <- c("inverse_psm_result", "psm_result", "list")
  result
}
