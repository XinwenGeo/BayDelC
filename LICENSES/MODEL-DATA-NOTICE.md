# Calibrated model and data notice

The files under `inst/models/main/` are reduced prediction assets derived from
the final BayDelC v0.3 calibrations. They contain posterior parameter draws or
posterior mean curves, calibration metadata, and the calibration table needed
for provenance. They do not contain Stan executables or full fitted `brms`
objects and cannot be used to recreate the original MCMC chains.

Users should cite the BayDelC article when using these calibrated assets.
Individual source observations retain the attribution recorded in each model's
`training_data` table and in the accompanying article. The software and bundled
assets are distributed under the package license; this does not alter rights in
the original publications from which calibration observations were compiled.
