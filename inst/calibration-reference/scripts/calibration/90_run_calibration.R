# BayDelC v0.3 clean calibration pipeline ----------------------------------
# Safe defaults do NOT refit or overwrite trained RDS files.
# To reproduce model training, explicitly change the relevant flag(s) to TRUE.

RUN_PREPARE_DATA <- TRUE
RUN_FIT_MAIN_MODELS <- FALSE
RUN_FIT_SENSITIVITY_MODELS <- FALSE
RUN_MODEL_EVALUATION <- TRUE
RUN_PUBLICATION_TABLES <- TRUE

if (RUN_PREPARE_DATA) {
  source("scripts/calibration/01_prepare_data.R")
}

if (RUN_FIT_MAIN_MODELS) {
  source("scripts/calibration/02_fit_main_models.R")
}

if (RUN_FIT_SENSITIVITY_MODELS) {
  source("scripts/calibration/03_fit_sensitivity_models.R")
}

if (RUN_MODEL_EVALUATION) {
  source("scripts/calibration/04_model_evaluation.R")
}

if (RUN_PUBLICATION_TABLES) {
  source("scripts/calibration/05_make_publication_tables.R")
}

message("BayDelC calibration pipeline completed according to RUN_* flags.")
