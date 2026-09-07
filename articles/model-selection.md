# Model and calibration selection

BLR estimates a monotonic linear relationship and is the preferred starting
point for transparent forward modelling and direct inverse uncertainty
propagation. GPR estimates a nonlinear relationship and is useful when that
flexibility is scientifically justified.

Inverse GPR results depend on both the likelihood and the prior. Always report
the prior, grid, observation-error representation, model/calibration choice,
and the posterior mass outside the calibration range. Do not select `Sub` or
`All` on predictive performance alone: the fossil taxa and analytical
definition must match the calibration.

## Primary and sensitivity collections

The published primary calibration is selected with `model_set = "main"`.
Sensitivity models are alternatives for robustness analysis rather than a pool
from which to select the most favorable answer:

- BLR `student_t` changes the residual distribution;
- GPR `shrinkage` adds the semiparametric shrinkage specification;
- GPR `latent_input` represents oxygen uncertainty through latent inputs and is
  available for the restricted-taxon calibration;
- GPR `x_jitter` preserves the 50-member diagnostic ensemble but is not a
  prediction-ready posterior model.

Use `baydelc_models("sensitivity", prediction_only = TRUE)` for valid modelling
choices. Report the collection and variant alongside every result.
