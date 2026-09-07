# Bundled model assets

These files are compact derivatives of the final BayDelC v0.3 model bundles.
They support prediction and inversion but do not contain fitted Stan objects,
raw MCMC chain files, log-likelihood matrices, or posterior predictive arrays.

- BLR assets contain `alpha`, `beta`, and residual-scale posterior draws.
- GPR assets contain 2,000 matched posterior mean curves evaluated from 0 to
  400 µmol kg⁻¹ at 1 µmol kg⁻¹ spacing, plus residual-scale draws.

`MODEL_MANIFEST.csv` records file size and MD5 checksum. The reproducible
extraction script is retained in `data-raw/build_portable_models.R` in the
source repository. Running that script does not refit or recalibrate a model.
