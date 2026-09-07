# BayDelC v0.3: fit the four main calibration models ------------------------
# The model functions, priors, seeds, datasets, and output paths are unchanged
# from the original scripts/02_fit_main_models.R.

source("scripts/00_setup.R")
source("config/model_registry.R")
source("R/data_prepare.R")
source("R/fit_blr.R")
source("R/fit_gpr.R")

prepared <- prepare_calibration_data(
  "data-raw/O2DD13C_NEW.xlsx",
  "data/processed"
)

sub_blr <- fit_blr_eiv_gaussian(
  prepared$calib_sub,
  "Sub_Linear_EIV_Gaussian",
  seed = 101
)

all_blr <- fit_blr_eiv_gaussian(
  prepared$calib_all,
  "All_Linear_EIV_Gaussian",
  seed = 102
)

sub_gpr_base <- fit_gpr_baseline(
  prepared$calib_sub,
  "Sub_GPR_Baseline",
  seed = 201,
  use_lscale_prior = TRUE
)

all_gpr_base <- fit_gpr_baseline(
  prepared$calib_all,
  "All_GPR_Baseline",
  seed = 202,
  use_lscale_prior = TRUE
)

save_model_bundle(
  sub_blr,
  baydelc_model_path("BLR_Sub_main", must_exist = FALSE)
)
save_model_bundle(
  all_blr,
  baydelc_model_path("BLR_All_main", must_exist = FALSE)
)
save_model_bundle(
  sub_gpr_base,
  baydelc_model_path("GPR_Sub_main", must_exist = FALSE)
)
save_model_bundle(
  all_gpr_base,
  baydelc_model_path("GPR_All_main", must_exist = FALSE)
)

message("Main models saved to inst/models/main")
