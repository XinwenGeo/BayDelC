# BayDelC v0.3: fit sensitivity calibration models --------------------------
# Defaults reproduce the original scripts/02b_fit_sensitivity_models.R.
# Set individual RUN_* flags to FALSE when only selected models are needed.

source("scripts/00_setup.R")
source("config/model_registry.R")
source("R/data_prepare.R")
source("R/fit_blr.R")
source("R/fit_gpr.R")
source("R/forward.R")

prepared <- prepare_calibration_data(
  "data-raw/O2DD13C_NEW.xlsx",
  "data/processed"
)

RUN_STUDENT_T <- TRUE
RUN_SHRINKAGE_GPR <- TRUE
RUN_XJITTER <- TRUE
RUN_LATENT_INPUT <- TRUE

if (RUN_STUDENT_T) {
  sub_blr_t <- fit_blr_eiv_student_t(
    prepared$calib_sub,
    "Sub_Linear_EIV_StudentT",
    seed = 111
  )
  all_blr_t <- fit_blr_eiv_student_t(
    prepared$calib_all,
    "All_Linear_EIV_StudentT",
    seed = 112
  )

  save_model_bundle(
    sub_blr_t,
    baydelc_model_path("BLR_Sub_StudentT", must_exist = FALSE)
  )
  save_model_bundle(
    all_blr_t,
    baydelc_model_path("BLR_All_StudentT", must_exist = FALSE)
  )
}

if (RUN_SHRINKAGE_GPR) {
  sub_gpr_shrink <- fit_gpr_shrinkage(
    prepared$calib_sub,
    "Sub_GPR_Shrinkage",
    seed = 211,
    use_lscale_prior = TRUE,
    b_positive = FALSE
  )
  all_gpr_shrink <- fit_gpr_shrinkage(
    prepared$calib_all,
    "All_GPR_Shrinkage",
    seed = 212,
    use_lscale_prior = TRUE,
    b_positive = FALSE
  )

  save_model_bundle(
    sub_gpr_shrink,
    baydelc_model_path("GPR_Sub_shrinkage", must_exist = FALSE)
  )
  save_model_bundle(
    all_gpr_shrink,
    baydelc_model_path("GPR_All_shrinkage", must_exist = FALSE)
  )
}

if (RUN_XJITTER) {
  # Original manuscript setting: M = 50, iter = 2000, warmup = 1000.
  sub_xjit <- fit_gpr_xjitter_ensemble(
    prepared$calib_sub,
    "Sub_GPR_XJitter_Baseline",
    M = 50,
    gpr_model = "baseline",
    seed = 3101,
    iter = 2000,
    warmup = 1000
  )
  all_xjit <- fit_gpr_xjitter_ensemble(
    prepared$calib_all,
    "All_GPR_XJitter_Baseline",
    M = 50,
    gpr_model = "baseline",
    seed = 3102,
    iter = 2000,
    warmup = 1000
  )

  saveRDS(
    sub_xjit,
    baydelc_model_path("GPR_Sub_xjitter", must_exist = FALSE)
  )
  saveRDS(
    all_xjit,
    baydelc_model_path("GPR_All_xjitter", must_exist = FALSE)
  )
}

if (RUN_LATENT_INPUT) {
  # Computationally intensive; original v0.3 fits the Sub model only.
  sub_latent <- fit_gpr_latent_input(
    prepared$calib_sub,
    "Sub_GPR_LatentInput",
    seed = 5101
  )

  save_model_bundle(
    sub_latent,
    baydelc_model_path("GPR_Sub_latent", must_exist = FALSE)
  )
}

message("Sensitivity models finished according to RUN_* flags.")
