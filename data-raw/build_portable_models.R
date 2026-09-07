# Build the small prediction assets shipped with BayDelC.
#
# This script NEVER refits a model. It reads the archived v0.3 fit bundles,
# extracts BLR parameter draws, and evaluates existing GPR posterior draws on a
# fixed BWO grid. Set BAYDELC_SOURCE_ROOT to override the adjacent source path.

source_root <- Sys.getenv(
  "BAYDELC_SOURCE_ROOT",
  unset = file.path("..", "BayDelC_v0.3_publishable")
)
source_dir <- file.path(normalizePath(source_root), "inst", "models", "main")
output_root <- file.path("inst", "models")
dir.create(file.path(output_root, "main"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "sensitivity"), recursive = TRUE, showWarnings = FALSE)

gpr_draw_count <- 2000L
gpr_grid <- seq(0, 400, by = 1)

even_draw_ids <- function(n, requested) {
  unique(as.integer(round(seq(1, n, length.out = min(n, requested)))))
}

portable_metadata <- function(x, source_file, backend) {
  list(
    model_name = x$model_name,
    model_type = x$model_type,
    algorithm = if (grepl("GPR", x$model_type, ignore.case = TRUE)) "GPR" else "BLR",
    bundle_version = x$bundle_version,
    data_version = x$data_version,
    training_data = x$training_data,
    data_filter_rule = x$data_filter_rule,
    calibration_range = x$calibration_range,
    settings = x$settings,
    notes = x$notes,
    source_created_time = x$created_time,
    prediction_backend = backend,
    portable_asset_version = "1",
    portable_created_with = paste(R.version$major, R.version$minor, sep = "."),
    source_file = basename(source_file)
  )
}

build_latent_gpr <- function(input, output) {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Building the portable latent-input GPR requires rstan, but using it does not.")
  }
  x <- readRDS(input)
  ext <- rstan::extract(x$fit_object)
  needed <- c("alpha", "eta", "rho", "sigma_model", "x_true", "f_latent")
  missing <- setdiff(needed, names(ext))
  if (length(missing)) stop("Latent GPR source lacks: ", paste(missing, collapse = ", "))
  ids <- even_draw_ids(length(ext$alpha), gpr_draw_count)
  tr <- x$settings$x_transform
  x_std <- (gpr_grid - as.numeric(tr$mean)) / as.numeric(tr$sd)
  curves <- matrix(NA_real_, nrow = length(ids), ncol = length(gpr_grid))
  for (k in seq_along(ids)) {
    i <- ids[[k]]
    xt <- as.numeric(ext$x_true[i, ])
    f <- as.numeric(ext$f_latent[i, ])
    eta <- as.numeric(ext$eta[i])
    rho <- as.numeric(ext$rho[i])
    covariance <- outer(xt, xt, function(a, b) eta^2 * exp(-0.5 * ((a - b) / rho)^2))
    diag(covariance) <- diag(covariance) + 1e-6
    cross_covariance <- outer(x_std, xt, function(a, b) eta^2 * exp(-0.5 * ((a - b) / rho)^2))
    curves[k, ] <- as.numeric(ext$alpha[i] + cross_covariance %*% solve(covariance, f))
  }
  out <- portable_metadata(x, input, "interpolated_latent_input_mean_curves")
  out$posterior_draws <- data.frame(sigma = abs(as.numeric(ext$sigma_model[ids])))
  out$portable_bwo_grid <- gpr_grid
  out$portable_mean_draws <- curves
  out$source_draw_ids <- ids
  saveRDS(out, output, compress = "xz", version = 3)
}

build_xjitter_diagnostic <- function(input, output) {
  x <- readRDS(input)
  x$algorithm <- "GPR"
  x$prediction_backend <- "diagnostic_curve_ensemble"
  x$prediction_ready <- FALSE
  saveRDS(x, output, compress = "xz", version = 3)
}

build_blr <- function(input, output) {
  x <- readRDS(input)
  keep <- intersect(c("alpha", "beta", "sigma_model", "nu"),
                    names(as.data.frame(x$posterior_draws)))
  if (!all(c("alpha", "beta") %in% keep)) {
    stop("BLR source lacks alpha or beta: ", input)
  }
  out <- portable_metadata(x, input, "posterior_parameter_draws")
  out$posterior_draws <- as.data.frame(x$posterior_draws)[keep]
  saveRDS(out, output, compress = "xz", version = 3)
}

build_gpr <- function(input, output) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Building portable GPR assets requires brms, but using them does not.")
  }
  x <- readRDS(input)
  all_draws <- as.data.frame(x$posterior_draws)
  ids <- even_draw_ids(nrow(all_draws), gpr_draw_count)
  tr <- x$settings$x_transform
  newdata <- data.frame(
    O2 = gpr_grid,
    O2_std = (gpr_grid - as.numeric(tr$mean)) / as.numeric(tr$sd),
    dd13C_sd = 0
  )
  # brms GP predictions conditionally sample function values at new points.
  # Fix this extraction seed so rebuilding the compact asset is byte-stable.
  extraction_seed <- if (is.null(x$settings$seed)) 90001L else as.integer(x$settings$seed) + 90000L
  set.seed(extraction_seed)
  curves <- as.matrix(brms::posterior_epred(
    x$fit_object, newdata = newdata, draw_ids = ids, re_formula = NA
  ))
  if (!all(dim(curves) == c(length(ids), length(gpr_grid)))) {
    stop("Unexpected GPR prediction dimensions for ", input)
  }
  sigma_name <- intersect(c("sigma", "sigma_model", "residual_sd"), names(all_draws))
  if (!length(sigma_name)) stop("GPR source lacks residual sigma: ", input)
  out <- portable_metadata(x, input, "interpolated_posterior_mean_curves")
  out$posterior_draws <- all_draws[ids, sigma_name[[1L]], drop = FALSE]
  out$portable_bwo_grid <- gpr_grid
  out$portable_mean_draws <- curves
  out$source_draw_ids <- ids
  saveRDS(out, output, compress = "xz", version = 3)
}

jobs <- list(
  list("BLR", "main", "Sub_Linear_EIV_Gaussian_v0.3.rds", "Sub_Linear_EIV_Gaussian_v0.3_portable.rds"),
  list("BLR", "main", "All_Linear_EIV_Gaussian_v0.3.rds", "All_Linear_EIV_Gaussian_v0.3_portable.rds"),
  list("GPR", "main", "Sub_GPR_Baseline_v0.3.rds", "Sub_GPR_Baseline_v0.3_portable.rds"),
  list("GPR", "main", "All_GPR_Baseline_v0.3.rds", "All_GPR_Baseline_v0.3_portable.rds"),
  list("BLR", "sensitivity", "Sub_Linear_EIV_StudentT_v0.3.rds", "Sub_Linear_EIV_StudentT_v0.3_portable.rds"),
  list("BLR", "sensitivity", "All_Linear_EIV_StudentT_v0.3.rds", "All_Linear_EIV_StudentT_v0.3_portable.rds"),
  list("GPR", "sensitivity", "Sub_GPR_Shrinkage_v0.3.rds", "Sub_GPR_Shrinkage_v0.3_portable.rds"),
  list("GPR", "sensitivity", "All_GPR_Shrinkage_v0.3.rds", "All_GPR_Shrinkage_v0.3_portable.rds"),
  list("XJITTER", "sensitivity", "Sub_GPR_XJitter_Baseline_v0.3.rds", "Sub_GPR_XJitter_Baseline_v0.3_diagnostic.rds"),
  list("XJITTER", "sensitivity", "All_GPR_XJitter_Baseline_v0.3.rds", "All_GPR_XJitter_Baseline_v0.3_diagnostic.rds"),
  list("LATENT", "sensitivity", "Sub_GPR_LatentInput_v0.3.rds", "Sub_GPR_LatentInput_v0.3_portable.rds")
)

for (job in jobs) {
  input <- file.path(normalizePath(source_root), "inst", "models", job[[2L]], job[[3L]])
  output <- file.path(output_root, job[[2L]], job[[4L]])
  message("Building ", basename(output), " without refitting")
  if (job[[1L]] == "BLR") build_blr(input, output)
  else if (job[[1L]] == "GPR") build_gpr(input, output)
  else if (job[[1L]] == "LATENT") build_latent_gpr(input, output)
  else build_xjitter_diagnostic(input, output)
}

files <- list.files(output_root, pattern = "[.]rds$", full.names = TRUE, recursive = TRUE)
manifest <- data.frame(
  file = substring(files, nchar(output_root) + 2L),
  bytes = unname(file.info(files)$size),
  md5 = unname(tools::md5sum(files)),
  stringsAsFactors = FALSE
)
utils::write.csv(manifest, file.path(output_root, "MODEL_MANIFEST.csv"), row.names = FALSE)
