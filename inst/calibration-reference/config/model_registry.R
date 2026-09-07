# BayDelC v0.3 active model registry ---------------------------------------
# This file contains paths and stable labels only. It does not load or fit a
# model unless baydelc_load_model() is called explicitly.

baydelc_model_registry <- function(root = ".") {
  rel_path <- c(
    "inst/models/main/Sub_Linear_EIV_Gaussian_v0.3.rds",
    "inst/models/main/All_Linear_EIV_Gaussian_v0.3.rds",
    "inst/models/main/Sub_GPR_Baseline_v0.3.rds",
    "inst/models/main/All_GPR_Baseline_v0.3.rds",
    "inst/models/sensitivity/Sub_Linear_EIV_StudentT_v0.3.rds",
    "inst/models/sensitivity/All_Linear_EIV_StudentT_v0.3.rds",
    "inst/models/sensitivity/Sub_GPR_Shrinkage_v0.3.rds",
    "inst/models/sensitivity/All_GPR_Shrinkage_v0.3.rds",
    "inst/models/sensitivity/Sub_GPR_LatentInput_v0.3.rds",
    "inst/models/sensitivity/Sub_GPR_XJitter_Baseline_v0.3.rds",
    "inst/models/sensitivity/All_GPR_XJitter_Baseline_v0.3.rds"
  )

  data.frame(
    model_id = c(
      "BLR_Sub_main", "BLR_All_main", "GPR_Sub_main", "GPR_All_main",
      "BLR_Sub_StudentT", "BLR_All_StudentT",
      "GPR_Sub_shrinkage", "GPR_All_shrinkage", "GPR_Sub_latent",
      "GPR_Sub_xjitter", "GPR_All_xjitter"
    ),
    model_name = c(
      "Sub_Linear_EIV_Gaussian", "All_Linear_EIV_Gaussian",
      "Sub_GPR_Baseline", "All_GPR_Baseline",
      "Sub_Linear_EIV_StudentT", "All_Linear_EIV_StudentT",
      "Sub_GPR_Shrinkage", "All_GPR_Shrinkage", "Sub_GPR_LatentInput",
      "Sub_GPR_XJitter_Baseline", "All_GPR_XJitter_Baseline"
    ),
    display_name = c(
      "Linear (Sub)", "Linear (All)", "GPR (Sub)", "GPR (All)",
      "Linear-t (Sub)", "Linear-t (All)",
      "GPR-shrinkage (Sub)", "GPR-shrinkage (All)", "GPR-latent (Sub)",
      "GPR-XJitter (Sub)", "GPR-XJitter (All)"
    ),
    short_name = c(
      "BLR (Sub)", "BLR (All)", "GPR (Sub)", "GPR (All)",
      "BLR-t (Sub)", "BLR-t (All)",
      "GPR-S (Sub)", "GPR-S (All)", "GPR-L (Sub)",
      "GPR-XJ (Sub)", "GPR-XJ (All)"
    ),
    algorithm = c(
      "BLR", "BLR", "GPR", "GPR", "BLR", "BLR",
      "GPR", "GPR", "GPR", "GPR", "GPR"
    ),
    calibration = c(
      "Sub", "All", "Sub", "All", "Sub", "All",
      "Sub", "All", "Sub", "Sub", "All"
    ),
    role = c(
      "main", "main", "main", "main",
      "sensitivity", "sensitivity", "sensitivity", "sensitivity",
      "sensitivity", "sensitivity", "sensitivity"
    ),
    backend = c(
      "BLR_EIV", "BLR_EIV", "brms_GPR", "brms_GPR",
      "BLR_EIV_StudentT", "BLR_EIV_StudentT",
      "brms_GPR", "brms_GPR", "rstan_latent_GPR",
      "GPR_ensemble", "GPR_ensemble"
    ),
    active = rep(TRUE, 11L),
    relative_path = rel_path,
    path = file.path(root, rel_path),
    stringsAsFactors = FALSE
  )
}

baydelc_model_path <- function(model_id, root = ".", must_exist = TRUE) {
  registry <- baydelc_model_registry(root = root)
  hit <- registry[registry$model_id == model_id, , drop = FALSE]
  if (nrow(hit) != 1L) {
    stop(sprintf("Unknown or duplicated model_id: %s", model_id), call. = FALSE)
  }
  path <- hit$path[[1L]]
  if (must_exist && !file.exists(path)) {
    stop(sprintf("Registered model file does not exist: %s", path), call. = FALSE)
  }
  path
}

baydelc_load_model <- function(model_id, root = ".") {
  readRDS(baydelc_model_path(model_id, root = root, must_exist = TRUE))
}

baydelc_main_model_ids <- function() {
  c("BLR_Sub_main", "BLR_All_main", "GPR_Sub_main", "GPR_All_main")
}

baydelc_sensitivity_model_ids <- function(include_xjitter = FALSE) {
  ids <- c(
    "BLR_Sub_StudentT", "BLR_All_StudentT",
    "GPR_Sub_shrinkage", "GPR_All_shrinkage", "GPR_Sub_latent"
  )
  if (isTRUE(include_xjitter)) {
    ids <- c(ids, "GPR_Sub_xjitter", "GPR_All_xjitter")
  }
  ids
}

baydelc_load_models <- function(model_ids, root = ".") {
  out <- lapply(model_ids, baydelc_load_model, root = root)
  names(out) <- model_ids
  out
}
