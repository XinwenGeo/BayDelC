
# BayDelC application data readers ------------------------------------------

find_age_col <- function(df) {
  nm <- names(df)
  hit <- nm[grep("^Age|Age.*ka|ka BP|age", nm, ignore.case = TRUE)][1]
  if (is.na(hit)) stop("No age column found. Columns are: ", paste(nm, collapse = ", "))
  hit
}

read_application_data <- function(data_ex = "DATA_EX.xlsx",
                                  ctrace_pmip = "ctrace-PMIP数据表.xlsx") {
  list(
    summary = readxl::read_excel(data_ex, sheet = "summary") |> dplyr::as_tibble(),
    data_ex_path = data_ex,
    ctrace = readxl::read_excel(ctrace_pmip, sheet = "ctrace") |> dplyr::as_tibble(),
    pmip = readxl::read_excel(ctrace_pmip, sheet = "PMIP") |> dplyr::as_tibble(),
    ctrace_pmip_path = ctrace_pmip
  )
}

get_site_metadata <- function(app, site) {
  app$summary |>
    dplyr::filter(.data$Site == site) |>
    dplyr::slice(1)
}

get_modern_bwo <- function(app, site) {
  x <- get_site_metadata(app, site)
  as.numeric(x[["Modern BWO"]][1])
}

get_site_depth <- function(app, site) {
  x <- get_site_metadata(app, site)
  as.numeric(x[["Data_Depth"]][1])
}

get_site_series <- function(app,
                            site,
                            max_age = NULL,
                            dd13c_sd_default = 0.05) {
  df <- readxl::read_excel(app$data_ex_path, sheet = site) |> dplyr::as_tibble()
  age_col <- find_age_col(df)
  if (!"dd13c" %in% names(df)) {
    ddcol <- names(df)[tolower(names(df)) %in% c("dd13c", "dd13c.", "dd13c_permil")][1]
    if (is.na(ddcol)) stop("No dd13c column found in sheet: ", site)
    names(df)[names(df) == ddcol] <- "dd13c"
  }
  out <- df |>
    dplyr::mutate(
      Age = as.numeric(.data[[age_col]]),
      dd13c = as.numeric(.data$dd13c),
      BWO_original = if ("BWO" %in% names(df)) as.numeric(.data$BWO) else NA_real_,
      dd13c_sd = dd13c_sd_default,
      Site = site
    ) |>
    dplyr::filter(is.finite(.data$Age), is.finite(.data$dd13c)) |>
    dplyr::arrange(.data$Age)

  if (!is.null(max_age)) out <- out |> dplyr::filter(.data$Age <= max_age)
  out
}

prepare_fig9_pmip_input <- function(app) {
  app$pmip |>
    dplyr::mutate(
      site_id = dplyr::row_number(),
      Site = as.character(.data$Site),
      Lon = as.numeric(.data$Data_long),
      Lat = as.numeric(.data$Data_Lat),
      Depth = as.numeric(.data$Data_Depth),
      O2 = as.numeric(.data$pmip_o2_mean_umolkg_approx),
      O2_sd = as.numeric(.data$pmip_o2_sd_umolkg_approx),
      obs_dd13c = as.numeric(.data$ts_mean_19_23ka),
      obs_dd13c_sd = as.numeric(.data$ts_sd_19_23ka)
    ) |>
    dplyr::filter(is.finite(.data$O2), is.finite(.data$obs_dd13c))
}

get_ctrace_series <- function(app, site, max_age = 21) {
  if (!site %in% names(app$ctrace)) {
    warning("Site not found in ctrace table: ", site)
    return(tibble::tibble(Age = numeric(), O2 = numeric(), Site = character()))
  }
  app$ctrace |>
    dplyr::transmute(Age = as.numeric(.data[["Age (ka)"]]),
                     O2 = as.numeric(.data[[site]]),
                     Site = site) |>
    dplyr::filter(is.finite(.data$Age), is.finite(.data$O2), .data$Age <= max_age) |>
    dplyr::arrange(.data$Age)
}

make_modern_marker <- function(app, site) {
  tibble::tibble(
    Site = site,
    Age = 0,
    O2 = get_modern_bwo(app, site),
    Depth = get_site_depth(app, site)
  )
}

interp_to_grid <- function(df, age_col = "Age", value_col = "O2", grid = seq(0, 21, by = 0.2)) {
  df <- df |> dplyr::filter(is.finite(.data[[age_col]]), is.finite(.data[[value_col]])) |> dplyr::arrange(.data[[age_col]])
  if (nrow(df) < 2) return(tibble::tibble(Age = grid, value = NA_real_))
  tibble::tibble(
    Age = grid,
    value = stats::approx(df[[age_col]], df[[value_col]], xout = grid, rule = 2, ties = "ordered")$y
  )
}

detrend_series <- function(df, age_col = "Age", value_col = "value",
                           span = 0.35) {
  d <- df |> dplyr::filter(is.finite(.data[[age_col]]), is.finite(.data[[value_col]]))
  if (nrow(d) < 6) {
    d$trend <- mean(d[[value_col]], na.rm = TRUE)
  } else {
    fit <- stats::loess(stats::as.formula(paste(value_col, "~", age_col)), data = d, span = span, degree = 2)
    d$trend <- stats::predict(fit, newdata = d)
  }
  d$residual <- d[[value_col]] - d$trend
  d
}

make_power_spectrum <- function(df, age_col = "Age", value_col = "value", dt = 0.2) {
  d <- df |> dplyr::filter(is.finite(.data[[age_col]]), is.finite(.data[[value_col]])) |> dplyr::arrange(.data[[age_col]])
  if (nrow(d) < 8) return(tibble::tibble(period = numeric(), power = numeric()))
  x <- scale(d[[value_col]])[, 1]
  sp <- stats::spec.pgram(x, taper = 0.1, log = "no", plot = FALSE, fast = FALSE)
  freq <- sp$freq / dt  # cycles per ka
  ok <- freq > 0
  period <- 1 / freq[ok]
  power <- sp$spec[ok]
  tibble::tibble(period = period, power = power / max(power, na.rm = TRUE))
}
