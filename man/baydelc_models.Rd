\name{baydelc_models}
\alias{baydelc_models}
\alias{baydelc_load_model}
\alias{baydelc_model_info}
\title{Inspect Bundled BayDelC Models}
\usage{
baydelc_models(model_set = c("all", "main", "sensitivity"),
  prediction_only = FALSE)

baydelc_load_model(algorithm = c("BLR", "GPR"),
  calibration = c("Sub", "All"),
  model_set = c("main", "sensitivity"), variant = NULL)

baydelc_model_info(algorithm = c("BLR", "GPR"),
  calibration = c("Sub", "All"),
  model_set = c("main", "sensitivity"), variant = NULL)
}
\arguments{
\item{algorithm}{Model family, either \code{"BLR"} or \code{"GPR"}.}
\item{calibration}{Calibration dataset, either restricted-taxon \code{"Sub"}
or mixed-taxon \code{"All"}.}
\item{model_set}{For listing, \code{"all"}, \code{"main"}, or
\code{"sensitivity"}; for loading, primary \code{"main"} or
\code{"sensitivity"}.}
\item{prediction_only}{Whether to omit diagnostic-only assets.}
\item{variant}{Optional model variant. Defaults depend on model collection and algorithm.}
}
\value{\code{baydelc_models()} returns a model registry data frame;
\code{baydelc_load_model()} returns a lightweight model bundle; and
\code{baydelc_model_info()} returns a named metadata list.}
\description{Lists, loads, or reports metadata for calibrated models shipped with BayDelC.}
