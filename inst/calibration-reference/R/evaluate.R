# Model evaluation and tables ----------------------------------------------

model_summary_table <- function(bundles) {
  purrr::map_dfr(bundles, function(b) {
    tibble::tibble(
      Model = b$model_name,
      Type = b$model_type,
      N = if (!is.null(b$training_data)) nrow(b$training_data) else NA_integer_,
      RMSE = b$metrics$RMSE %||% NA_real_,
      Bayesian_R2 = b$metrics$Bayesian_R2 %||% NA_real_,
      ELPD = b$metrics$ELPD %||% extract_elpd(b$loo_object, "Estimate"),
      ELPD_SE = b$metrics$ELPD_SE %||% extract_elpd(b$loo_object, "SE"),
      max_Pareto_k = b$metrics$max_pareto_k %||% extract_pareto_k(b$loo_object),
      PPC_coverage_95 = b$metrics$PPC_coverage_95 %||% NA_real_,
      Calibration_min_O2 = min(b$training_data$O2, na.rm = TRUE),
      Calibration_max_O2 = max(b$training_data$O2, na.rm = TRUE),
      Bundle_version = b$bundle_version %||% NA_character_
    )
  })
}

write_model_tables <- function(main_paths, sensitivity_paths = character(), outdir = "outputs/tables") {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  main_bundles <- load_model_bundles(main_paths)
  main_summary <- model_summary_table(main_bundles)
  readr::write_csv(main_summary, file.path(outdir, "Table_main_model_comparison_v0.3.csv"))
  readr::write_csv(posterior_summary_table(main_bundles), file.path(outdir, "Table_main_posterior_summary_v0.3.csv"))

  if (length(sensitivity_paths) > 0) {
    sens_bundles <- load_model_bundles(sensitivity_paths)
    sens_summary <- model_summary_table(sens_bundles)
    readr::write_csv(sens_summary, file.path(outdir, "TableS_model_comparison_sensitivity_v0.3.csv"))
    readr::write_csv(posterior_summary_table(sens_bundles), file.path(outdir, "TableS_posterior_summary_sensitivity_v0.3.csv"))
  }
  invisible(TRUE)
}
