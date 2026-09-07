# BayDelC v0.3: evaluate existing model bundles -----------------------------
# This script never fits a model. It reads the registered RDS files and writes
# the same technical comparison tables as the original v0.3 evaluation script.

source("scripts/00_setup.R")
source("config/model_registry.R")
source("R/evaluate.R")

main_paths <- vapply(
  baydelc_main_model_ids(),
  baydelc_model_path,
  character(1),
  must_exist = TRUE
)

sensitivity_ids <- baydelc_sensitivity_model_ids(include_xjitter = TRUE)
sensitivity_paths <- vapply(
  sensitivity_ids,
  baydelc_model_path,
  character(1),
  must_exist = FALSE
)
sensitivity_paths <- sensitivity_paths[file.exists(sensitivity_paths)]

# X-jitter ensembles do not have the scalar bundle structure used by
# write_model_tables(); they are summarized separately below.
sensitivity_paths_scalar <- sensitivity_paths[
  !grepl("XJitter", sensitivity_paths, ignore.case = TRUE)
]

write_model_tables(
  main_paths,
  sensitivity_paths_scalar,
  outdir = "outputs/tables"
)

xjitter_paths <- sensitivity_paths[
  grepl("XJitter", sensitivity_paths, ignore.case = TRUE)
]

if (length(xjitter_paths) > 0L) {
  xjit_tables <- purrr::map_dfr(xjitter_paths, function(path) {
    bundle <- readRDS(path)
    bundle$member_metrics |>
      dplyr::mutate(
        Model = bundle$model_name,
        Type = bundle$model_type
      )
  })

  readr::write_csv(
    xjit_tables,
    "outputs/tables/TableS_xjitter_member_metrics_v0.3.csv"
  )
}

message("Technical model tables saved to outputs/tables")
