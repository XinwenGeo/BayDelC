# Sensitivity model assets

- `*_StudentT_*`: BLR calibrations with Student-t residuals.
- `*_Shrinkage_*`: semiparametric GPR sensitivity calibrations.
- `*_LatentInput_*`: GPR calibration treating oxygen inputs as latent.
- `*_XJitter_*`: 50-member oxygen-jitter ensemble curves for diagnostics only.

Use `baydelc_models("sensitivity")` to see supported combinations and whether an
asset is prediction-ready.
