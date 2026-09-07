# BayDelC v0.3: prepare the fixed calibration datasets ----------------------
# Scientific behavior is unchanged from scripts/01_prepare_data.R.

source("scripts/00_setup.R")
source("R/data_prepare.R")

prepared <- prepare_calibration_data(
  "data-raw/O2DD13C_NEW.xlsx",
  "data/processed"
)

readr::write_csv(
  prepared$decision_log,
  "outputs/tables/TableS_data_decision_log_v0.3.csv"
)

message(
  "Data prepared: ALL n=", nrow(prepared$calib_all),
  "; SUB n=", nrow(prepared$calib_sub)
)
