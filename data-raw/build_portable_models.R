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
output_dir <- file.path("inst", "models", "main")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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
  list("BLR", "Sub_Linear_EIV_Gaussian_v0.3.rds", "Sub_Linear_EIV_Gaussian_v0.3_portable.rds"),
  list("BLR", "All_Linear_EIV_Gaussian_v0.3.rds", "All_Linear_EIV_Gaussian_v0.3_portable.rds"),
  list("GPR", "Sub_GPR_Baseline_v0.3.rds", "Sub_GPR_Baseline_v0.3_portable.rds"),
  list("GPR", "All_GPR_Baseline_v0.3.rds", "All_GPR_Baseline_v0.3_portable.rds")
)

for (job in jobs) {
  input <- file.path(source_dir, job[[2L]])
  output <- file.path(output_dir, job[[3L]])
  message("Building ", basename(output), " without refitting")
  if (job[[1L]] == "BLR") build_blr(input, output) else build_gpr(input, output)
}

files <- list.files(output_dir, pattern = "[.]rds$", full.names = TRUE)
manifest <- data.frame(
  file = basename(files),
  bytes = unname(file.info(files)$size),
  md5 = unname(tools::md5sum(files)),
  stringsAsFactors = FALSE
)
utils::write.csv(manifest, file.path(output_dir, "MODEL_MANIFEST.csv"), row.names = FALSE)
