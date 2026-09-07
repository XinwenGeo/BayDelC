#' Prepare versioned BayDelC calibration data
#' @param xlsx Path to O2DD13C_NEW.xlsx.
#' @param outdir Output directory.
#' @return A list with raw_data, calib_all, calib_sub, and decision_log.
prepare_calibration_data <- function(xlsx = "data-raw/O2DD13C_NEW.xlsx", outdir = "data/processed") {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  raw <- readxl::read_excel(xlsx, sheet = "Raw_Data")
  dat <- raw |>
    dplyr::mutate(
      O2 = as.numeric(.data[["Average oxygen (umol/kg)"]]),
      O2_sd = as.numeric(.data[["Average oxygen (umol/kg) 1stdev"]]),
      dd13C = as.numeric(.data[["Ddl3C"]]),
      dd13C_sd = as.numeric(.data[["Ddl3C 1stdev"]]),
      Lon = as.numeric(.data[["Lon"]]),
      Lat = as.numeric(.data[["Lat"]]),
      water_depth = as.numeric(.data[["Water Depth"]]),
      NO3 = suppressWarnings(as.numeric(.data[["[NO3][µmol/kg]"]])),
      PO4 = suppressWarnings(as.numeric(.data[["[PO4][µmol/kg]"]])),
      outlier = dplyr::coalesce(as.numeric(.data[["outlier"]]), 0),
      Group = dplyr::case_when(
        grepl("bottom and pore water", .data$type, ignore.case = TRUE) ~ "BPW",
        grepl("Globobulimina|G\\.affinis|G\\.pacifica", .data$type, ignore.case = TRUE) & !grepl("epifaunal", .data$type, ignore.case = TRUE) ~ "BWG",
        grepl("epifaunal", .data$type, ignore.case = TRUE) & .data$`Cib Species` == "C.wuellerstorfi" ~ "C_wuel",
        grepl("epifaunal", .data$type, ignore.case = TRUE) & .data$`Cib Species` == "C.pachyderma" ~ "C_pach",
        grepl("epifaunal", .data$type, ignore.case = TRUE) & .data$`Cib Species` == "C.mundulus" ~ "C_mund",
        grepl("epifaunal", .data$type, ignore.case = TRUE) & grepl("unidentified", .data$`Cib Species`, ignore.case = TRUE) ~ "C_unid",
        TRUE ~ "Other"
      ),
      is_foram_only = grepl("epifaunal and deep infaunal benthic foraminifera", .data$type, ignore.case = TRUE),
      is_outlier = .data$outlier == 1,
      low_o2 = .data$O2 < 50,
      is_sub = .data$is_foram_only & .data$`Cib Species` == "C.wuellerstorfi",
      complete_core = !is.na(.data$O2) & !is.na(.data$O2_sd) & !is.na(.data$dd13C) & !is.na(.data$dd13C_sd),
      use_all_calib = .data$is_foram_only & !.data$is_outlier & !.data$low_o2 & .data$complete_core,
      use_sub_calib = .data$is_sub & !.data$is_outlier & !.data$low_o2 & .data$complete_core,
      calibration_decision = dplyr::case_when(
        .data$use_all_calib ~ "included_all",
        !.data$is_foram_only ~ "not_foram_only",
        .data$is_outlier ~ "outlier",
        .data$low_o2 ~ "O2_below_50",
        !.data$complete_core ~ "missing_core_variable",
        TRUE ~ "excluded_other"
      )
    )
  calib_all <- dplyr::filter(dat, .data$use_all_calib)
  calib_sub <- dplyr::filter(dat, .data$use_sub_calib)
  decision_log <- dplyr::select(dat, Location, `Core name`, Lat, Lon, O2, O2_sd, `Cib Species`, `Globobulimina species`, dd13C, dd13C_sd, type, reference, Group, is_outlier, low_o2, is_foram_only, use_all_calib, use_sub_calib, calibration_decision)
  readr::write_csv(dat, file.path(outdir, "raw_data_v0.2.csv"))
  readr::write_csv(calib_all, file.path(outdir, "calib_all_o2ge50_v0.2.csv"))
  readr::write_csv(calib_sub, file.path(outdir, "calib_sub_o2ge50_v0.2.csv"))
  readr::write_csv(decision_log, file.path(outdir, "data_decision_log_v0.2.csv"))
  list(raw_data = dat, calib_all = calib_all, calib_sub = calib_sub, decision_log = decision_log)
}
