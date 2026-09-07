# BayDelC v0.3: publication-facing model tables -----------------------------
# The table calculations are identical to the former publication plotting
# helper, but no plotting file or plotting patch is sourced here.

source("scripts/00_setup.R")
source("config/model_registry.R")
source("R/model_bundle.R")
source("R/evaluate.R")

model_nice_names <- c(
  Sub_Linear_EIV_Gaussian = "Linear (Sub)",
  All_Linear_EIV_Gaussian = "Linear (All)",
  Sub_GPR_Baseline = "GPR (Sub)",
  All_GPR_Baseline = "GPR (All)",
  Sub_Linear_EIV_StudentT = "Linear-t (Sub)",
  All_Linear_EIV_StudentT = "Linear-t (All)",
  Sub_GPR_Shrinkage = "GPR-shrinkage (Sub)",
  All_GPR_Shrinkage = "GPR-shrinkage (All)",
  Sub_GPR_LatentInput = "GPR-latent (Sub)"
)

nice_model_name <- function(x) {
  dplyr::recode(x, !!!model_nice_names, .default = x)
}

main_ids <- baydelc_main_model_ids()
sensitivity_ids <- baydelc_sensitivity_model_ids(include_xjitter = FALSE)

main_paths <- vapply(
  main_ids,
  baydelc_model_path,
  character(1),
  must_exist = TRUE
)

sensitivity_paths <- vapply(
  sensitivity_ids,
  baydelc_model_path,
  character(1),
  must_exist = FALSE
)
sensitivity_paths <- sensitivity_paths[file.exists(sensitivity_paths)]

main_bundles <- load_model_bundles(main_paths)
sensitivity_bundles <- if (length(sensitivity_paths) > 0L) {
  load_model_bundles(sensitivity_paths)
} else {
  list()
}

main_summary <- model_summary_table(main_bundles)
main_summary$Model_publication <- nice_model_name(main_summary$Model)

readr::write_csv(
  main_summary,
  "outputs/tables/Table1_main_model_comparison.csv"
)
readr::write_csv(
  posterior_summary_table(main_bundles),
  "outputs/tables/TableS1_main_posterior_summary.csv"
)

if (length(sensitivity_bundles) > 0L) {
  sensitivity_summary <- model_summary_table(sensitivity_bundles)
  sensitivity_summary$Model_publication <- nice_model_name(
    sensitivity_summary$Model
  )

  readr::write_csv(
    sensitivity_summary,
    "outputs/tables/TableS2_sensitivity_model_comparison.csv"
  )
  readr::write_csv(
    posterior_summary_table(sensitivity_bundles),
    "outputs/tables/TableS3_sensitivity_posterior_summary.csv"
  )
}

message("Publication tables saved to outputs/tables")
