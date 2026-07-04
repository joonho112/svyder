# Extract Draws Matrix from a svyder Object

Returns the corrected draws if available, otherwise the original draws.
This method allows `svyder` objects to be used wherever a numeric matrix
of posterior draws is expected.

## Usage

``` r
# S3 method for class 'svyder'
as.matrix(x, ...)
```

## Arguments

- x:

  A `svyder` object.

- ...:

  Ignored.

## Value

A numeric matrix of posterior draws (S x d).

## References

Lee, J., Williams, M. R., & Savitsky, T. D. (2026). Design Effect Ratios
for Bayesian Survey Models: A Diagnostic Framework for Identifying
Survey-Sensitive Parameters. *Journal of Survey Statistics and
Methodology*. Submitted.

## See also

[`der_correct()`](https://joonho112.github.io/svyder/reference/der_correct.md)
for applying the correction.

Other svyder-methods:
[`glance.svyder()`](https://joonho112.github.io/svyder/reference/glance.svyder.md),
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
  cluster = nsece_demo$psu, family = "binomial",
  sigma_theta = nsece_demo$sigma_theta,
  param_types = nsece_demo$param_types
)
draws <- as.matrix(result)
dim(draws)
#> [1] 4000   54
```
