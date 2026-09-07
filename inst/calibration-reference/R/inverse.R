# BayDelC inverse model -----------------------------------------------------
# Inverse maps observed Δδ13C to posterior BWO/O2 distributions.

normalize_weights <- function(w) {
  s <- sum(w, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) rep(1 / length(w), length(w)) else w / s
}

weighted_grid_quantile <- function(grid, weights, probs = c(0.025, 0.5, 0.975)) {
  weights <- normalize_weights(weights)
  cdf <- cumsum(weights)
  stats::setNames(vapply(probs, function(p) grid[which.min(abs(cdf - p))], numeric(1)), probs)
}

inverse_blr <- function(bundle, dd13C_obs, dd13C_sd = 0,
                        prior_mean = NULL, prior_sd = NULL,
                        prior_lower = 0, prior_upper = 350,
                        mode = c("strict", "practical", "deterministic"),
                        O2_grid = seq(0, 350, by = 0.5),
                        n_draw = 2000, seed = 9101,
                        return_density = FALSE) {
  mode <- match.arg(mode)
  stopifnot(grepl("BLR_EIV", bundle$model_type))
  draws <- as.data.frame(bundle$posterior_draws)
  draws <- sample_draws_df(draws, n_draw, seed)

  if (mode == "practical") {
    draws <- data.frame(
      alpha = stats::median(draws$alpha, na.rm = TRUE),
      beta = stats::median(draws$beta, na.rm = TRUE),
      sigma_model = stats::median(draws$sigma_model, na.rm = TRUE)
    )
  }
  if (mode == "deterministic") {
    alpha <- stats::median(draws$alpha, na.rm = TRUE)
    beta <- stats::median(draws$beta, na.rm = TRUE)
    med <- (dd13C_obs - alpha) / beta
    return(tibble::tibble(sample = seq_along(dd13C_obs), median = med, q025 = NA_real_, q975 = NA_real_, mode = mode))
  }

  grid <- O2_grid[O2_grid >= prior_lower & O2_grid <= prior_upper]
  base_prior <- rep(1, length(grid))
  if (!is.null(prior_mean) && !is.null(prior_sd)) {
    prior_mean <- rep(prior_mean, length.out = length(dd13C_obs))
    prior_sd <- rep(prior_sd, length.out = length(dd13C_obs))
  }
  dd13C_sd <- rep(dd13C_sd, length.out = length(dd13C_obs))

  all_density <- list()
  res <- lapply(seq_along(dd13C_obs), function(i) {
    prior_i <- base_prior
    if (!is.null(prior_mean) && !is.null(prior_sd)) {
      prior_i <- stats::dnorm(grid, prior_mean[i], prior_sd[i])
    }
    w_draw <- matrix(NA_real_, nrow = nrow(draws), ncol = length(grid))
    for (j in seq_len(nrow(draws))) {
      mu <- draws$alpha[j] + draws$beta[j] * grid
      sd_j <- sqrt(draws$sigma_model[j]^2 + dd13C_sd[i]^2)
      w_draw[j, ] <- normalize_weights(stats::dnorm(dd13C_obs[i], mu, sd_j) * prior_i)
    }
    wbar <- normalize_weights(colMeans(w_draw, na.rm = TRUE))
    qs <- weighted_grid_quantile(grid, wbar)
    if (return_density) all_density[[i]] <<- tibble::tibble(sample = i, O2 = grid, density = wbar)
    tibble::tibble(sample = i, median = qs["0.5"], q025 = qs["0.025"], q975 = qs["0.975"], mode = mode)
  }) |>
    dplyr::bind_rows()

  if (return_density) attr(res, "density") <- dplyr::bind_rows(all_density)
  res
}

inverse_bwo <- function(bundle, ...) inverse_blr(bundle, ...)
