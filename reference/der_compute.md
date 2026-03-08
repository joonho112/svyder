# Compute Design Effect Ratios

Computes Design Effect Ratios (DER) for each parameter in a Bayesian
hierarchical model fitted to complex survey data. The DER measures how
much each parameter's posterior variance should change when the survey
design is properly accounted for.

## Usage

``` r
der_compute(x, ...)

# S3 method for class 'matrix'
der_compute(
  x,
  ...,
  y,
  X,
  group,
  weights,
  psu = NULL,
  family = "binomial",
  sigma_theta,
  sigma_e = NULL,
  beta_prior_sd = 5,
  param_types = NULL,
  design = NULL
)

# Default S3 method
der_compute(x, ...)

# S3 method for class 'brmsfit'
der_compute(
  x,
  ...,
  weights,
  group = NULL,
  psu = NULL,
  family = NULL,
  sigma_theta = NULL,
  sigma_e = NULL,
  beta_prior_sd = 5,
  param_types = NULL,
  design = NULL
)

# S3 method for class 'CmdStanMCMC'
der_compute(
  x,
  ...,
  y,
  X,
  group,
  weights,
  psu = NULL,
  family = "binomial",
  sigma_theta,
  sigma_e = NULL,
  beta_prior_sd = 5,
  param_types = NULL,
  design = NULL
)

# S3 method for class 'stanreg'
der_compute(
  x,
  ...,
  weights,
  psu = NULL,
  sigma_theta = NULL,
  sigma_e = NULL,
  beta_prior_sd = 5,
  param_types = NULL,
  design = NULL
)
```

## Arguments

- x:

  A draws matrix (S x d), brmsfit, CmdStanMCMC, or stanreg object. For
  the matrix method, columns 1:p are fixed effects and columns
  (p+1):(p+J) are random effects.

- ...:

  Additional arguments passed to methods.

- y:

  Response vector (length N).

- X:

  Design matrix (N x p).

- group:

  Integer group indicator (1 to J).

- weights:

  Survey weights (positive, length N).

- psu:

  PSU indicators (default: same as group).

- family:

  Model family: `"binomial"` or `"gaussian"`.

- sigma_theta:

  Estimated random effect SD.

- sigma_e:

  Residual SD (gaussian only).

- beta_prior_sd:

  Prior SD for fixed effects (default: 5).

- param_types:

  Character vector of length p: `"fe_within"` or `"fe_between"`.

- design:

  A `survey.design2` object (alternative to weights/psu).

## Value

A `svyder` object containing DER values, sandwich matrices, per-group
diagnostics, and original posterior draws.

## Details

DER is defined as the ratio of the sandwich variance (which accounts for
survey design) to the naive MCMC posterior variance. Values near 1
indicate that the posterior is already well-calibrated; values
substantially above 1 indicate underestimation of uncertainty.

**brmsfit method**: The brms method auto-detects the response vector
`y`, design matrix `X`, grouping variable `group`, model `family`, and
random effect SD `sigma_theta` from the fitted model object. The user
must provide `weights` (survey weights) and optionally `psu` (primary
sampling unit indicators), since these are not stored in brms fit
objects.

The `group` argument can optionally be provided to override
auto-detection from the random effects structure.

**CmdStanMCMC method**: CmdStan does not store model data in the fit
object, so the user must provide all data arguments: `y`, `X`, `group`,
`weights`, `sigma_theta`, and optionally `psu`, `sigma_e` (for
gaussian), `beta_prior_sd`, and `param_types`.

The draws are extracted from the CmdStanMCMC object using
[`extract_draws()`](https://joonho112.github.io/svyder/reference/extract_draws.md)
and then passed to the matrix method. The draws matrix must have columns
ordered as `[beta_1, ..., beta_p, theta_1, ..., theta_J]`.

**stanreg method**: The rstanarm method auto-detects the response vector
`y`, design matrix `X`, grouping variable `group`, model `family`,
random effect SD `sigma_theta`, and residual SD `sigma_e` (for gaussian)
from the fitted model object. The user must provide `weights` (survey
weights) and optionally `psu` (primary sampling unit indicators).

## See also

[`der_classify()`](https://joonho112.github.io/svyder/reference/der_classify.md)
for classification,
[`der_correct()`](https://joonho112.github.io/svyder/reference/der_correct.md)
for correction,
[`der_diagnose()`](https://joonho112.github.io/svyder/reference/der_diagnose.md)
for the all-in-one pipeline.

Other core-pipeline:
[`der_classify()`](https://joonho112.github.io/svyder/reference/der_classify.md),
[`der_correct()`](https://joonho112.github.io/svyder/reference/der_correct.md),
[`der_diagnose()`](https://joonho112.github.io/svyder/reference/der_diagnose.md)

## Examples

``` r
data(nsece_demo)
result <- der_compute(
  nsece_demo$draws,
  y = nsece_demo$y, X = nsece_demo$X,
  group = nsece_demo$group, weights = nsece_demo$weights,
  psu = nsece_demo$psu, family = "binomial",
  sigma_theta = nsece_demo$sigma_theta,
  param_types = nsece_demo$param_types
)
print(result)
#> svyder diagnostic (54 parameters)
#>   Family: binomial | N = 6785 | J = 51
#>   DER range: [0.235, 5.315]
#>   (not yet classified -- run der_classify())
#>   Compute time: 0.188 sec
```
