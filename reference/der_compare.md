# Compare DER Across Clustering Definitions

Computes DER using different clustering (PSU) definitions and compares
the results. This is useful for comparing, e.g., state-level vs
PSU-level clustering to assess sensitivity of DER to the choice of
primary sampling unit.

## Usage

``` r
der_compare(x, clusters, ...)
```

## Arguments

- x:

  A draws matrix (S x d) or a `svyder` object. If a `svyder` object, the
  original draws and data are extracted automatically.

- clusters:

  A named list of PSU vectors to compare. Each element should be an
  integer vector of length N.

- ...:

  Additional arguments passed to
  [`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md).
  Required when `x` is a matrix: `y`, `X`, `group`, `weights`, `family`,
  `sigma_theta`, etc.

## Value

A `data.frame` with columns: `param`, `cluster_name`, `der`.

## See also

[`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md)
for computing DER with a single clustering.

Other analysis:
[`der_decompose()`](https://joonho112.github.io/svyder/reference/der_decompose.md),
[`der_sensitivity()`](https://joonho112.github.io/svyder/reference/der_sensitivity.md),
[`der_theorem_check()`](https://joonho112.github.io/svyder/reference/der_theorem_check.md)

## Examples

``` r
data(nsece_demo)
# Compare DER using original PSU vs group-level clustering
comp <- der_compare(
  nsece_demo$draws,
  clusters = list(
    psu   = nsece_demo$psu,
    group = nsece_demo$group
  ),
  y = nsece_demo$y, X = nsece_demo$X,
  group = nsece_demo$group, weights = nsece_demo$weights,
  family = "binomial",
  sigma_theta = nsece_demo$sigma_theta,
  param_types = nsece_demo$param_types
)
head(comp)
#>      param cluster_name       der
#> 1  beta[1]          psu 0.2617696
#> 2  beta[2]          psu 2.6868825
#> 3  beta[3]          psu 0.3427294
#> 4 theta[1]          psu 3.3838218
#> 5 theta[2]          psu 0.6758803
#> 6 theta[3]          psu 1.1187639
```
