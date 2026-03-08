# Glance at a svyder Object

Returns a one-row data frame with model-level summary statistics
following broom conventions. Provides an overview of the DER diagnostic
results.

## Usage

``` r
glance.svyder(x, ...)
```

## Arguments

- x:

  A `svyder` object.

- ...:

  Additional arguments (unused).

## Value

A `data.frame` with one row and columns:

- n_params:

  Total number of parameters.

- n_flagged:

  Number of parameters flagged for correction.

- pct_flagged:

  Percentage of parameters flagged.

- tau:

  Classification threshold.

- family:

  Model family (`"binomial"` or `"gaussian"`).

- n_obs:

  Number of observations.

- n_groups:

  Number of groups/clusters.

- mean_deff:

  Mean per-group design effect.

- mean_B:

  Mean per-group shrinkage factor.

- der_min:

  Minimum DER value.

- der_max:

  Maximum DER value.

## See also

[`tidy.svyder()`](https://joonho112.github.io/svyder/reference/tidy.svyder.md)
for per-parameter summaries.

Other svyder-methods:
[`as.matrix.svyder()`](https://joonho112.github.io/svyder/reference/as.matrix.svyder.md),
[`is.svyder()`](https://joonho112.github.io/svyder/reference/is.svyder.md),
[`print.svyder()`](https://joonho112.github.io/svyder/reference/print.svyder.md),
[`summary.svyder()`](https://joonho112.github.io/svyder/reference/summary.svyder.md),
[`tidy.svyder()`](https://joonho112.github.io/svyder/reference/tidy.svyder.md)

## Examples

``` r
data(nsece_demo)
result <- der_diagnose(
  nsece_demo$draws,
  y = nsece_demo$y, X = nsece_demo$X,
  group = nsece_demo$group, weights = nsece_demo$weights,
  psu = nsece_demo$psu, family = "binomial",
  sigma_theta = nsece_demo$sigma_theta,
  param_types = nsece_demo$param_types
)
glance.svyder(result)
#>   n_params n_flagged pct_flagged tau   family n_obs n_groups mean_deff
#> 1       54        30    55.55556 1.2 binomial  6785       51   2.59527
#>      mean_B   der_min  der_max
#> 1 0.8543053 0.2349819 5.314837
```
