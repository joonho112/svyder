# Compute Design Effect Ratios

Computes Design Effect Ratios (DER) for each parameter in a Bayesian
hierarchical model fitted to complex survey data. The DER of a parameter
is the ratio of its design-based sandwich variance – computed for a
*declared variance target* – to its model-based MCMC posterior variance.
Values substantially above 1 indicate that the posterior understates
design-based uncertainty; values at or below 1 indicate that correction
is unnecessary (and typically harmful).

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
  cluster,
  strata = NULL,
  cluster_strata = NULL,
  family = "binomial",
  sigma_theta,
  sigma_e = NULL,
  normalize = c("unit_mean", "group_size", "none"),
  center_meat = TRUE,
  df_correct = TRUE,
  beta_prior_sd = Inf,
  param_types = NULL,
  design = NULL,
  psu = NULL
)

# Default S3 method
der_compute(x, ...)

# S3 method for class 'brmsfit'
der_compute(
  x,
  ...,
  weights,
  group = NULL,
  cluster = NULL,
  strata = NULL,
  psu = NULL,
  family = NULL,
  sigma_theta = NULL,
  sigma_e = NULL,
  normalize = c("unit_mean", "group_size", "none"),
  center_meat = TRUE,
  df_correct = TRUE,
  beta_prior_sd = Inf,
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
  cluster = NULL,
  strata = NULL,
  psu = NULL,
  family = "binomial",
  sigma_theta,
  sigma_e = NULL,
  normalize = c("unit_mean", "group_size", "none"),
  center_meat = TRUE,
  df_correct = TRUE,
  beta_prior_sd = Inf,
  param_types = NULL,
  design = NULL
)

# S3 method for class 'stanreg'
der_compute(
  x,
  ...,
  weights,
  cluster = NULL,
  strata = NULL,
  psu = NULL,
  sigma_theta = NULL,
  sigma_e = NULL,
  normalize = c("unit_mean", "group_size", "none"),
  center_meat = TRUE,
  df_correct = TRUE,
  beta_prior_sd = Inf,
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

- cluster:

  Sandwich aggregation unit (length N). **Required** (unless supplied
  via `design`): pass the design PSU ids for a design-based target or
  `group` for a model-aligned target.

- strata:

  Optional stratum indicator (length N).

- cluster_strata:

  Optional vector giving the stratum of *every* cluster in the declared
  design universe, indexed by cluster id (length G, the universe size).
  Supplying it declares the cluster universe explicitly, so cluster ids
  with no sampled observations – selected PSUs whose second-stage (e.g.,
  Poisson) sample is empty – enter the meat as zero score-total
  clusters, contributing to the stratum centering and the
  \\C_h/(C_h-1)\\ correction rather than silently vanishing from the
  target. When supplied, `cluster` must contain integer ids in
  `1:length(cluster_strata)` and `strata` is required. Leave `NULL` when
  every cluster of the realized design has at least one observation
  (e.g., the NSECE application).

- family:

  Model family: `"binomial"` or `"gaussian"`.

- sigma_theta:

  Estimated random effect SD (plug-in).

- sigma_e:

  Residual SD (gaussian only).

- normalize:

  Weight convention: `"unit_mean"` (default), `"group_size"`, or
  `"none"`. Recorded in the result.

- center_meat:

  Center cluster score totals within strata (default `TRUE`). Set
  `FALSE` only to reproduce v1 results.

- df_correct:

  Apply the \\C_h/(C_h-1)\\ stratum correction (default `TRUE`). Set
  `FALSE` only to reproduce v1 results.

- beta_prior_sd:

  Prior SD for fixed effects in the bread (default `Inf` = diffuse-beta
  convention, no ridge).

- param_types:

  Character vector of length p: `"fe_within"` or `"fe_between"`. If
  `NULL`, types are inferred from the design matrix (columns constant
  within every group are `"fe_between"`) and the inference is reported
  via a message.

- design:

  A `survey.design2` object; supplies `weights`, `cluster`
  (first-stage), and `strata` (first-stage).

- psu:

  Deprecated alias for `cluster`.

## Value

A `svyder` object (a structured list) with, among others:

- der:

  Named numeric vector of Design Effect Ratios, one per parameter
  (`beta[1..p]` then `theta[1..J]`).

- V_sand:

  Sandwich (declared-target) variance matrix \\H^{-1} J_c H^{-1}\\ (\\d
  \times d\\).

- sigma_mcmc:

  Model-based MCMC posterior covariance (\\d \times d\\); `der` is its
  diagonal ratio to `V_sand`.

- target:

  List recording the declared variance target: number of clusters,
  strata, centering / DF-correction flags, weight convention, and bread
  convention.

- excluded:

  Data frame of Tier III parameters (hyperparameters) excluded from the
  DER diagnostic, with reasons.

- data:

  List of data slots (`X`, `group`, `weights`, `cluster`, `strata`,
  working weights `v`, score residuals `r`) retained for re-targeting
  and decomposition.

- classification:

  Data frame with one row per parameter (`param_name`, `param_type`,
  `der`); tiers are filled in by
  [`der_classify()`](https://joonho112.github.io/svyder/reference/der_classify.md).

- deff_j, B_j:

  Per-group Kish design effects and shrinkage factors (length J).

## Details

**brmsfit method**: The brms method auto-detects the response vector
`y`, design matrix `X`, grouping variable `group`, model `family`, and
random effect SD `sigma_theta` from the fitted model object. The user
must provide `weights` (survey weights) and `cluster` (the sandwich
aggregation unit), since these are not stored in brms fit objects.
`strata`, `normalize`, `center_meat`, and `df_correct` are passed
through to the matrix method.

The `group` argument can optionally be provided to override
auto-detection from the random effects structure.

**CmdStanMCMC method**: CmdStan does not store model data in the fit
object, so the user must provide all data arguments: `y`, `X`, `group`,
`weights`, `cluster`, `sigma_theta`, and optionally `strata`, `sigma_e`
(for gaussian), `normalize`, `center_meat`, `df_correct`,
`beta_prior_sd`, and `param_types`.

The draws are extracted from the CmdStanMCMC object using
[`extract_draws()`](https://joonho112.github.io/svyder/reference/extract_draws.md)
and then passed to the matrix method. The draws matrix must have columns
ordered as `[beta_1, ..., beta_p, theta_1, ..., theta_J]`.

**stanreg method**: The rstanarm method auto-detects the response vector
`y`, design matrix `X`, grouping variable `group`, model `family`,
random effect SD `sigma_theta`, and residual SD `sigma_e` (for gaussian)
from the fitted model object. The user must provide `weights` (survey
weights) and `cluster` (the sandwich aggregation unit); `strata`,
`normalize`, `center_meat`, and `df_correct` are passed through to the
matrix method.

## The declared variance target

The sandwich variance is not unique: it depends on choices that are part
of the estimand, not implementation details. `der_compute()` therefore
requires them explicitly and records them in the returned object:

- `cluster` – the aggregation unit of the score outer product (the
  design PSU for a design-based target, or the model group for a
  model-aligned target). There is no default: DER values can change
  materially between the two (in the NSECE application the flagged set
  moves from 1 of 54 to 28 of 54, including 27 of the 51 state effects).

- `strata` – optional design strata; cluster totals are centered within
  strata and a \\C_h/(C_h-1)\\ correction is applied.

- `normalize` – the survey-weight convention. Random-effect DERs are not
  invariant to weight rescaling, so the convention is part of the
  target. `"unit_mean"` (default) rescales weights to mean 1 overall;
  `"group_size"` rescales within each group to sum to the group size;
  `"none"` uses the weights as supplied.

## The bread

The bread is the penalized (posterior) Hessian at the posterior mean:
the pseudo-likelihood curvature plus the random-effect prior precision
\\1/\sigma\_\theta^2\\. The likelihood-only Hessian is singular for any
model with an intercept and a full set of group effects, so the prior
term is what makes the sandwich well-defined. By the diffuse-beta
convention the fixed-effect prior contributes no curvature by default
(`beta_prior_sd = Inf`); supply a finite value to reproduce pipelines
that used an informative beta prior in the bread.

## References

Lee, J., Williams, M. R., & Savitsky, T. D. (2026). Design Effect Ratios
for Bayesian Survey Models: A Diagnostic Framework for Identifying
Survey-Sensitive Parameters. *Journal of Survey Statistics and
Methodology*. Submitted.

Savitsky, T. D., & Williams, M. R. (2022). Pseudo Bayesian mixed models
under informative sampling. *Journal of Official Statistics*, 38(4),
901–928.
[doi:10.2478/jos-2022-0039](https://doi.org/10.2478/jos-2022-0039)

Williams, M. R., & Savitsky, T. D. (2021). Uncertainty estimation for
pseudo-Bayesian inference under complex sampling. *International
Statistical Review*, 89(1), 72–107.
[doi:10.1111/insr.12376](https://doi.org/10.1111/insr.12376)

Binder, D. A. (1983). On the variances of asymptotically normal
estimators from complex surveys. *International Statistical Review*,
51(3), 279–292. [doi:10.2307/1402588](https://doi.org/10.2307/1402588)

Kish, L. (1965). *Survey Sampling*. John Wiley & Sons.

## See also

[`der_classify()`](https://joonho112.github.io/svyder/reference/der_classify.md)
for classification,
[`der_correct()`](https://joonho112.github.io/svyder/reference/der_correct.md)
for correction,
[`der_compare()`](https://joonho112.github.io/svyder/reference/der_compare.md)
for comparing targets,
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
  cluster = nsece_demo$psu, family = "binomial",
  sigma_theta = nsece_demo$sigma_theta,
  param_types = nsece_demo$param_types
)
print(result)
#> svyder diagnostic (54 parameters)
#>   Family: binomial | N = 6785 | J = 51
#>   Target: 644 cluster(s) | meat: centered, DF-corrected | weights: unit_mean
#>   DER range: [0.235, 5.308]
#>   Tier III (DER undefined): sigma_theta
#>   (not yet classified -- run der_classify())
#>   Compute time: 0.019 sec
```
