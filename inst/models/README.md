# Bundled model assets

`main/` contains the four primary v0.3 calibrations. `sensitivity/` contains
Student-t BLR, shrinkage GPR, latent-input GPR, and X-jitter curve diagnostics.
`MODEL_MANIFEST.csv` records the size and checksum of every asset.

Prediction-ready BLR files retain posterior parameter draws. Prediction-ready
GPR files retain 2,000 posterior mean curves on a 0–400 µmol kg⁻¹ grid and
matched residual-scale draws. X-jitter assets retain the original 50-member
curve summaries and are diagnostic-only. No asset contains a compiled Stan
object, full fitted model, log-likelihood matrix, or raw MCMC chain.
