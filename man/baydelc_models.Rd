\name{baydelc_models}
\alias{baydelc_models}
\alias{baydelc_load_model}
\alias{baydelc_model_info}
\title{Inspect Bundled BayDelC Models}
\usage{
baydelc_models()

baydelc_load_model(algorithm = c("BLR", "GPR"),
  calibration = c("Sub", "All"))

baydelc_model_info(algorithm = c("BLR", "GPR"),
  calibration = c("Sub", "All"))
}
\arguments{
\item{algorithm}{Model family, either \code{"BLR"} or \code{"GPR"}.}
\item{calibration}{Calibration dataset, either restricted-taxon \code{"Sub"}
or mixed-taxon \code{"All"}.}
}
\value{\code{baydelc_models()} returns a model registry data frame;
\code{baydelc_load_model()} returns a lightweight model bundle; and
\code{baydelc_model_info()} returns a named metadata list.}
\description{Lists, loads, or reports metadata for calibrated models shipped with BayDelC.}
