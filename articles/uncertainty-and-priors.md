# Uncertainty and priors

BayDelC keeps distinct uncertainty sources explicit. Forward calculations can
propagate uncertain bottom-water oxygen (BWO) inputs, calibration uncertainty,
and residual proxy-model variability. Inverse calculations propagate
calibration uncertainty, a selected representation of Δδ¹³C observation
uncertainty, residual variability, and, when requested, an external BWO prior.

All uncertainty arguments in this article are one-standard-deviation values
unless stated otherwise. BWO is expressed in µmol kg⁻¹ and Δδ¹³C in ‰.
BayDelC uses the epifaunal-minus-infaunal sign convention.

## Forward uncertainty

`forward_psm()` accepts fixed BWO values, Gaussian BWO uncertainty, empirical
BWO draws, or a custom input distribution. The default
`prediction = "predictive"` includes residual calibration scatter and therefore
describes a new proxy observation. Use `prediction = "latent"` when the target
is the uncertain calibrated mean relationship rather than a new observation.

```r
library(BayDelC)

# Uncertain BWO inputs and a predictive Δδ¹³C distribution.
fwd_predictive <- forward_psm(
  bwo = c(80, 150, 220),
  bwo_sd = c(5, 10, 10),
  algorithm = "BLR",
  calibration = "Sub",
  prediction = "predictive",
  n_draw = 2000,
  seed = 42,
  return_draws = TRUE
)

# The same inputs, excluding residual variability.
fwd_latent <- forward_psm(
  bwo = c(80, 150, 220),
  bwo_sd = c(5, 10, 10),
  algorithm = "BLR",
  calibration = "Sub",
  prediction = "latent",
  n_draw = 2000,
  seed = 42
)
```

The difference between these results is intentional: calibration-parameter
uncertainty is present in both, whereas residual variability is present only in
the predictive calculation.

## Inverse observation uncertainty

`inverse_psm()` accepts five mutually exclusive representations of Δδ¹³C
uncertainty. With `error_method = "auto"`, BayDelC selects the representation
implied by the supplied arguments. Specify `error_method` explicitly in an
archived analysis so that the intended calculation is unambiguous.

### Supplied total standard deviation

Use `dd13c_sd` when the complete observation-level uncertainty has already been
calculated. This is the simplest and usually the preferred interface.

```r
inv_total <- inverse_psm(
  dd13c = c(1.10, 1.65),
  dd13c_sd = c(0.08, 0.10),
  error_method = "provided_total",
  algorithm = "BLR",
  calibration = "Sub",
  n_draw = 2000,
  seed = 42,
  return_draws = TRUE
)
```

### H15-style default proxy error

When no study-specific uncertainty is available, the H15-style representation
uses

$$
s_{obs} = \sqrt{0.05^2 + 0.05^2 + 0.05^2}\ \text{‰}.
$$

It must be selected deliberately in reproducible analyses rather than treated
as a universal property of the proxy.

```r
inv_h15 <- inverse_psm(
  dd13c = 1.50,
  error_method = "h15_proxy_default",
  algorithm = "BLR",
  calibration = "Sub",
  n_draw = 2000,
  seed = 42
)
```

### Lu-style site variability

For the site-variability representation, BayDelC combines the supplied
infaunal site standard deviation with a 0.05‰ analytical component:

$$
s_{obs} = \sqrt{0.05^2 + s_{site}^2}.
$$

```r
inv_site <- inverse_psm(
  dd13c = c(1.25, 1.60),
  glob_site_sd = c(0.07, 0.11),
  error_method = "lu_site",
  algorithm = "BLR",
  calibration = "Sub",
  n_draw = 2000,
  seed = 42
)
```

### Component-wise uncertainty and covariance

If epifaunal and infaunal contributions are available separately, the total
variance is

$$
s_{obs}^2 = s_{epi}^2 + s_{infa}^2
             - 2\rho s_{epi}s_{infa} + s_{pairing}^2 + s_{other}^2.
$$

Each taxon-specific term can itself combine measurement and sampling
components in quadrature. The covariance term prevents correlated endmember
errors from being counted twice.

```r
err <- psm_error_components(
  epi_measurement_sd = 0.04,
  epi_sample_sd = 0.03,
  infa_measurement_sd = 0.04,
  infa_sample_sd = 0.05,
  pairing_sd = 0.02,
  rho_epi_infa = 0.25
)

inv_components <- inverse_psm(
  dd13c = 1.50,
  error_components = err,
  error_method = "components",
  algorithm = "BLR",
  calibration = "Sub",
  n_draw = 2000,
  seed = 42
)
```

Values supplied with `scale = "se"` are treated as standard deviations of the
corresponding observation-level estimates; BayDelC does not convert a raw
sample standard deviation to a standard error without a sample size.

### Empirical observation draws

Use empirical draws when observation uncertainty is asymmetric, non-Gaussian,
or already represented by an ensemble. Matrices always use observations in
rows and joint realizations in columns.

```r
set.seed(42)
dd13c_ensemble <- matrix(
  stats::rnorm(2000, mean = 1.50, sd = 0.08),
  nrow = 1
)

inv_draws <- inverse_psm(
  dd13c = 1.50,
  dd13c_draws = dd13c_ensemble,
  error_method = "draws",
  algorithm = "BLR",
  calibration = "Sub",
  n_draw = 2000,
  seed = 42,
  return_draws = TRUE
)
```

For several observations, a column should represent one joint realization
across all rows. Optional `dd13c_draw_weights` must have one value per column.

## Strict BLR inversion without a prior

For BLR with `prior = "none"`, the default `inversion = "auto"` selects strict
draw-wise inversion. For posterior draw \(s\), BayDelC evaluates

$$
x_i^{(s)} =
\frac{y_i^{(s)} - \epsilon_{model,i}^{(s)} - \alpha^{(s)}}
     {\beta^{(s)}}.
$$

Joint posterior rows preserve covariance among the calibration parameters.
With `include_model_residual = TRUE`, the residual draw is included, as in the
principal article applications. Strict inversion is not clipped to physical or
calibration bounds; implausible tails are retained and reported as diagnostics.

```r
inv_strict <- inverse_psm(
  dd13c = 1.50,
  dd13c_sd = 0.08,
  algorithm = "BLR",
  calibration = "Sub",
  inversion = "strict",
  prior = "none",
  include_model_residual = TRUE,
  n_draw = 4000,
  seed = 42,
  return_draws = TRUE
)
```

Setting `include_model_residual = FALSE` answers a different question by
excluding unexplained calibration scatter. Report this choice explicitly.

## Grid inversion and priors

BLR and GPR can both evaluate the inverse likelihood on a BWO grid. For a
candidate BWO value \(x\), observation \(y_i\), and calibration draw \(s\), the
Gaussian likelihood uses

$$
y_i \sim \mathcal{N}\left(\mu_s(x),
\sqrt{s_{obs,i}^2 + \sigma_{model,s}^2}\right).
$$

BayDelC averages this likelihood over calibration draws, multiplies by the BWO
prior and grid-cell width, and normalizes the resulting posterior mass. The
default grid spans `o2_bounds = c(0, 300)` at `o2_step = 1`.

### Physical uniform prior

This prior assigns constant density from 0 to 300 µmol kg⁻¹ and zero density
outside that interval.

```r
inv_gpr <- inverse_psm(
  dd13c = 1.50,
  dd13c_sd = 0.08,
  algorithm = "GPR",
  calibration = "Sub",
  inversion = "grid",
  prior = "uniform_physical",
  o2_bounds = c(0, 300),
  integration_draws = 500,
  n_draw = 2000,
  seed = 42,
  return_draws = TRUE,
  return_density = TRUE
)
```

### Normal or truncated-Normal prior

A Normal prior can represent independent external information. Supplying
`bounds` truncates the prior to the stated interval.

```r
inv_normal_prior <- inverse_psm(
  dd13c = 1.50,
  dd13c_sd = 0.08,
  algorithm = "BLR",
  calibration = "Sub",
  inversion = "grid",
  prior = "normal",
  prior_args = list(mean = 120, sd = 30, bounds = c(0, 300)),
  include_model_residual = TRUE,
  integration_draws = 500,
  n_draw = 2000,
  seed = 42,
  return_draws = TRUE
)
```

For a time series, `prior_args$mean` and `prior_args$sd` may contain one value
per observation, allowing an independently justified time-varying prior.

### Ensemble prior draws

An ensemble can preserve a non-Gaussian or jointly varying prior. As with
observation draws, rows are observations and columns are joint realizations.

```r
set.seed(42)
prior_ensemble <- matrix(
  pmin(300, pmax(0, stats::rnorm(3000, mean = 120, sd = 30))),
  nrow = 1
)

inv_ensemble_prior <- inverse_psm(
  dd13c = 1.50,
  dd13c_sd = 0.08,
  algorithm = "BLR",
  calibration = "Sub",
  inversion = "grid",
  prior = "draws",
  prior_args = list(draws = prior_ensemble),
  integration_draws = 500,
  n_draw = 2000,
  seed = 42,
  return_draws = TRUE
)
```

Numerical densities and custom density functions are also accepted through
`prior = "density"` and `prior = "custom"`. These should be normalized only
through BayDelC; the supplied values need to be non-negative but do not need to
integrate to one beforehand.

## Interpreting diagnostics

Inverse summaries report `P_below_zero`, `P_below_50`, `P_above_300`, and
`P_outside_calibration`. These probabilities expose weak constraint,
extrapolation, or sensitivity to physical limits; they are not automatic data
filters. A posterior concentrated at a grid boundary indicates that the grid or
prior bounds materially affect the result.

The reported `q025` and `q975` values are equal-tailed 95% uncertainty limits.
Their interpretation depends on the selected observation-error representation,
whether residual variability was included, the calibration model, and any BWO
prior.

## Reproducible reporting

Record at least the following with every archived analysis:

- package and model version;
- `algorithm`, `calibration`, `model_set`, and `variant`;
- `error_method` and the source and scale of its inputs;
- `inversion`, `prior`, `prior_args`, and grid bounds;
- `include_model_residual`, `n_draw`, `integration_draws`, and `seed`;
- posterior mass outside the calibration range and near imposed bounds.

A prior is part of the scientific model. Do not choose it after inspecting the
desired result. When more than one prior is scientifically plausible, report a
sensitivity analysis and retain a no-prior or weak-reference calculation where
the selected model permits it.
