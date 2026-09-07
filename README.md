# BayDelC

<!-- badges: start -->
[![R-CMD-check](https://github.com/clairezhang1mol/BayDelC/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/clairezhang1mol/BayDelC/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

BayDelC is an R package for forward and inverse proxy-system modelling of
bottom-water oxygen (BWO) using the benthic foraminiferal carbon-isotope offset
Δδ¹³C. It distributes the final BayDelC v0.3 calibrations as compact prediction
assets. Users can therefore run the PSM without installing Stan, rerunning MCMC,
or downloading the large raw chain files.

The two public modelling functions are:

- `forward_psm()`: BWO → Δδ¹³C;
- `inverse_psm()`: Δδ¹³C → a posterior distribution for BWO.

## Installation

Install the local source package from the directory containing
`BayDelC_package`:

```r
install.packages("BayDelC_package", repos = NULL, type = "source")
```

After the repository is published, users will be able to install it with:

```r
# install.packages("remotes")
remotes::install_github("clairezhang1mol/BayDelC")
```

BayDelC requires only base-recommended R packages for normal use. `brms` and
`posterior` are optional and are needed only when supplying certain custom full
model objects.

## Five-minute example

```r
library(BayDelC)

# See all bundled calibration models.
baydelc_models()

# Forward PSM: expected proxy values for three BWO inputs.
fwd <- forward_psm(
  bwo = c(80, 150, 220),
  bwo_sd = c(5, 10, 10),
  algorithm = "BLR",
  calibration = "Sub",
  n_draw = 2000,
  seed = 42,
  return_draws = TRUE
)
summary(fwd)
plot(fwd)

# Inverse PSM: reconstruct BWO from Δδ13C and its total 1-sigma error.
inv <- inverse_psm(
  dd13c = c(1.10, 1.65),
  dd13c_sd = c(0.08, 0.10),
  algorithm = "BLR",
  calibration = "Sub",
  n_draw = 2000,
  seed = 42,
  return_draws = TRUE
)
summary(inv)
plot(inv)
```

For nonlinear inversion, select GPR and use a physically explicit prior:

```r
gpr <- inverse_psm(
  dd13c = 1.50,
  dd13c_sd = 0.08,
  algorithm = "GPR",
  calibration = "Sub",
  inversion = "grid",
  prior = "uniform_physical",
  o2_bounds = c(0, 300),
  integration_draws = 500,
  n_draw = 2000,
  return_density = TRUE
)
```

## Choosing a model

BayDelC contains four primary models and seven sensitivity assets. List all
choices with:

```r
baydelc_models("all")
baydelc_models("sensitivity", prediction_only = TRUE)
```

The primary models are:

| Algorithm | `Sub` calibration | `All` calibration |
|---|---|---|
| BLR | Bayesian linear regression using the restricted taxon set | Bayesian linear regression using the mixed taxon set |
| GPR | Nonlinear Gaussian-process calibration using the restricted taxon set | Nonlinear Gaussian-process calibration using the mixed taxon set |

Use `Sub` when the target observations match the restricted taxonomic
definition used in that calibration. Use `All` only when the mixed-taxon
calibration matches the data. BLR is the transparent default and supports
direct inversion. GPR represents nonlinear structure but should be inverted
with an explicit BWO prior because nonlinear inverse mappings can be broad or
multimodal.

Inspect exact model identity and calibration bounds with:

```r
baydelc_model_info("BLR", "Sub")
```

Select a sensitivity model with `model_set` and `variant`:

```r
# Gaussian main BLR (the default)
main <- inverse_psm(
  dd13c = 1.5, dd13c_sd = 0.08,
  algorithm = "BLR", calibration = "Sub", model_set = "main"
)

# Student-t BLR sensitivity calibration
student_t <- inverse_psm(
  dd13c = 1.5, dd13c_sd = 0.08,
  algorithm = "BLR", calibration = "Sub",
  model_set = "sensitivity", variant = "student_t"
)

# GPR sensitivity calibrations
gpr_shrinkage <- forward_psm(
  bwo = c(80, 150, 220), algorithm = "GPR", calibration = "Sub",
  model_set = "sensitivity", variant = "shrinkage"
)
gpr_latent <- forward_psm(
  bwo = c(80, 150, 220), algorithm = "GPR", calibration = "Sub",
  model_set = "sensitivity", variant = "latent_input"
)
```

If `variant` is omitted, sensitivity BLR defaults to `student_t` and
sensitivity GPR defaults to `shrinkage`. The `x_jitter` entries preserve the
article's 50-member diagnostic curve ensembles. Because the archived X-jitter
objects do not contain posterior fits, they are intentionally marked
`prediction_ready = FALSE`; inspect them with `baydelc_load_model()` rather
than passing them to the forward or inverse PSM.

## Units and sign convention

- BWO is in µmol kg⁻¹.
- Δδ¹³C is in ‰.
- BayDelC v0.3 uses **epifaunal minus infaunal** δ¹³C. Reverse-signed input
  must be multiplied by −1 before use.
- Reported `q025` and `q975` are equal-tailed 95% uncertainty limits.

The package reports `P_outside_calibration`, `P_below_zero`, `P_above_300`, and
`P_below_50` to make extrapolation and physically sensitive posterior mass
visible. These are diagnostics, not automatic data filters. Scientific
interpretation should consider preservation, ecological offsets, chronology,
and whether the selected calibration is appropriate for the sample assemblage.

## Observation uncertainty

If a total Δδ¹³C standard deviation is available, pass it as `dd13c_sd`. If the
epifaunal and infaunal contributions are separate, use:

```r
err <- psm_error_components(
  epi_measurement_sd = 0.04,
  epi_sample_sd = 0.03,
  infa_measurement_sd = 0.04,
  infa_sample_sd = 0.05,
  rho_epi_infa = 0
)

inverse_psm(
  dd13c = 1.50,
  error_components = err,
  error_method = "components"
)
```

For non-Gaussian observation or prior uncertainty, supply an
observation-by-realization matrix through `dd13c_draws`, or prior ensemble
draws through `prior_args$draws`. Matrices always use observations in rows and
joint realizations in columns.

## Reproducibility and model assets

The bundled BLR assets retain posterior parameter draws. The bundled GPR assets
retain posterior mean functions evaluated on a dense BWO grid plus matched
residual-scale draws. This preserves posterior predictive use while removing
platform-specific compiled Stan objects. The package also includes calibration
source code, Stan programs, the small calibration workbook, and processed
tables for transparency. Full fitted objects and raw chains remain in the
associated article reproducibility archive.

Use `set.seed()` or the functions' `seed` argument for repeatable Monte Carlo
results. See `system.file("examples", package = "BayDelC")` for complete scripts
and `citation("BayDelC")` for the current citation record.

## Reproducing article figures as package examples

Two package-native examples reproduce the structure and scientific content of
article Figures 2 and 3 using only the distributed compact models:

```r
example_dir <- system.file("examples", "paper", package = "BayDelC")
source(file.path(example_dir, "Figure02_calibration_curves.R"))
source(file.path(example_dir, "Figure03_BLR_prediction.R"))
```

They write PDF and PNG files to `BayDelC_paper_examples/`, or to the directory
set in the `BAYDELC_EXAMPLE_OUTPUT` environment variable. They require the
optional packages `ggplot2` and `patchwork`. These are package-native
reproductions and may differ by small Monte Carlo or rendering details from the
submission files. The paper repository remains the authority for pixel-exact
publication output.

The complete main-versus-sensitivity usage example is available at:

```r
source(system.file("examples", "03_main_and_sensitivity.R", package = "BayDelC"))
```

## Calibration reference scripts

The original v0.3 calibration workflow is installed for expert inspection:

```r
calibration_reference <- system.file("calibration-reference", package = "BayDelC")
file.show(file.path(calibration_reference, "README.md"))
```

It includes data preparation, main and sensitivity fitting scripts, Stan model
files, the 28 KB source workbook, processed calibration tables, and the model
registry. Nothing in this directory runs when BayDelC is installed or loaded.
Recalibration is computationally expensive and creates a new model version; copy
the directory into a separate writable project and read its README before use.

## Scope

BayDelC is a calibrated scientific model, not a general oxygen sensor and not a
substitute for site-specific quality control. Predictions outside the
calibration domain are extrapolations. Full calibration provenance and
scientific limitations are described in the companion article.

## Development

Contributions are welcome through issues and focused pull requests. See
`CONTRIBUTING.md`. The package code and compact model assets are released under
the MIT license. Citation metadata intentionally marks the article as “in
preparation” until final journal details and DOI are available.
