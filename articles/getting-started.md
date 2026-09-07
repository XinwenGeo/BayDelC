# Getting started

## 1. Check the sign and units

Use BWO in µmol kg⁻¹ and Δδ¹³C in ‰. BayDelC expects Δδ¹³C = epifaunal δ¹³C −
infaunal δ¹³C.

## 2. Select the calibration

Choose `calibration = "Sub"` when the observations follow the restricted-taxon
definition; otherwise use `"All"` only when the mixed-taxon calibration is
scientifically appropriate. `baydelc_models()` lists the installed assets.

## 3. Run and retain diagnostics

```r
library(BayDelC)
x <- inverse_psm(
  dd13c = c(1.10, 1.50, 2.10),
  dd13c_sd = 0.08,
  algorithm = "BLR",
  calibration = "Sub",
  n_draw = 4000,
  seed = 2026,
  return_draws = TRUE
)
x$summary
x$diagnostics
```

Archive the package version (`baydelc_version()`), model metadata, arguments,
seed, summary, and diagnostics with each analysis.
