# BayDelC 0.3.0

- First public package release candidate.
- Adds unified `forward_psm()` and `inverse_psm()` interfaces.
- Bundles primary and sensitivity BLR posterior parameters and portable GPR
  prediction curves for the restricted-taxon (`Sub`) and mixed-taxon (`All`)
  calibrations.
- Adds explicit `model_set` and `variant` selection, including Student-t BLR,
  shrinkage GPR, latent-input GPR, and X-jitter diagnostic assets.
- Adds package-native Figure 2 and Figure 3 reproduction examples and the
  complete calibration workflow as an inert reference directory.
- Supports fixed, normal, draw-based, and custom uncertainty propagation.
- Adds diagnostic probabilities for physically or calibration-relevant ranges.
