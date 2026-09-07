# BayDelC calibration reference

This directory preserves the scripts, Stan programs, and calibration input
table used to fit the BayDelC v0.3 main and sensitivity models. It is included
for transparency and expert reuse. Installing or loading the BayDelC package
does **not** execute any of these files.

The fitted model assets distributed with the package are already calibrated.
Ordinary users should call `forward_psm()` or `inverse_psm()` and should not
refit the models.

## Contents

- `data-raw/O2DD13C_NEW.xlsx`: calibration workbook used by the workflow;
- `data/processed/`: resulting calibration tables and decision log;
- `R/`: data preparation, fitting, evaluation, forward and inverse helpers;
- `inst/stan/`: Gaussian BLR, Student-t BLR, and latent-input GPR programs;
- `scripts/calibration/`: preparation, main-model fitting, sensitivity fitting,
  evaluation, publication-table, and orchestration scripts.

## Inspecting the files after installation

```r
ref <- system.file("calibration-reference", package = "BayDelC")
list.files(ref, recursive = TRUE)
```

## Recalibration

Copy this directory to a separate writable project folder and use that copy as
the working directory. Install the packages listed in `REQUIREMENTS.txt`, then
review every seed, prior, sampler control, path, and `RUN_*` flag before running:

```r
source("scripts/calibration/90_run_calibration.R")
```

The safe defaults prepare data and evaluate existing full RDS files, but those
full RDS files are intentionally not included in this lightweight package.
Set `RUN_MODEL_EVALUATION <- FALSE` and `RUN_PUBLICATION_TABLES <- FALSE` when
starting from the supplied source inputs. Enabling either fitting flag performs
computationally intensive MCMC and writes new model bundles. It is a new
calibration run and will not alter the compact models installed with BayDelC.

The original v0.3 models were calibrated with the versions recorded inside
their bundles. Exact historical environment recreation belongs to the complete
article reproducibility archive; package updates should validate new fits
against the published version before replacing any distributed asset.
