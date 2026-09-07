# Uncertainty and priors

`forward_psm()` separates uncertain BWO inputs from residual proxy-model
variability. Set `prediction = "latent"` to describe only the calibrated mean
relationship, or keep `"predictive"` for a new proxy observation.

`inverse_psm()` can use total Gaussian measurement error, decomposed error
components, or empirical observation draws. For grid inversion, the available
priors include a physical uniform prior, normal or truncated-normal priors,
empirical draws, a supplied density, and a custom density function.

A prior is part of the scientific model. Avoid choosing it after inspecting the
desired result. Use sensitivity analyses when more than one prior is plausible,
and report posterior mass near grid limits as a warning that bounds matter.
