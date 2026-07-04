# Compare DER Across Variance Targets

Computes DER under different sandwich aggregation units (e.g., the
design PSU vs the model group) and compares the results side by side.
Because the aggregation unit is part of the variance-target definition,
reporting both is the recommended practice whenever the design PSUs and
the model groups differ.

## Usage

``` r
der_compare(x, clusters, ...)
```

## Arguments

- x:

  A `svyder` object (preferred) or a draws matrix (S x d).

- clusters:

  A named list of cluster-id vectors (each length N), one per target to
  compare.

- ...:

  For the matrix path only: arguments passed to
  [`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md)
  (`y`, `X`, `group`, `weights`, `family`, `sigma_theta`, ...).

## Value

A `data.frame` with columns: `param`, `cluster_name`, `der`.

## Details

When `x` is a `svyder` object from svyder \>= 0.2.0, the stored data
slots are reused: only the meat is rebuilt per target (the bread and the
MCMC covariance are target-invariant), so the comparison is fast and
requires no re-supply of `y`/`X`/`weights`. If the stored object used
`cluster_strata` to declare selected-but-empty clusters, that declared
universe is preserved when a comparison target is exactly the original
stored target.

## References

Lee, J., Williams, M. R., & Savitsky, T. D. (2026). Design Effect Ratios
for Bayesian Survey Models: A Diagnostic Framework for Identifying
Survey-Sensitive Parameters. *Journal of Survey Statistics and
Methodology*. Submitted.

## See also

[`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md)
for computing DER under a single target.

Other analysis:
[`der_decompose()`](https://joonho112.github.io/svyder/reference/der_decompose.md),
[`der_sensitivity()`](https://joonho112.github.io/svyder/reference/der_sensitivity.md),
[`der_theorem_check()`](https://joonho112.github.io/svyder/reference/der_theorem_check.md)

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
comp <- der_compare(result, clusters = list(
  design_psu  = nsece_demo$psu,
  model_group = nsece_demo$group
))
head(comp)
#>      param cluster_name       der
#> 1  beta[1]   design_psu 0.2626885
#> 2  beta[2]   design_psu 2.6894416
#> 3  beta[3]   design_psu 0.3430527
#> 4 theta[1]   design_psu 3.3858856
#> 5 theta[2]   design_psu 0.6759045
#> 6 theta[3]   design_psu 1.1182861
```
