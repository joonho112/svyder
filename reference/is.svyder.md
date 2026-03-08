# Test if an Object is a svyder Object

Checks whether an R object inherits from class `"svyder"`.

## Usage

``` r
is.svyder(x)
```

## Arguments

- x:

  An R object.

## Value

Logical; `TRUE` if `x` inherits from class `"svyder"`, `FALSE`
otherwise.

## See also

Other svyder-methods:
[`as.matrix.svyder()`](https://joonho112.github.io/svyder/reference/as.matrix.svyder.md),
[`glance.svyder()`](https://joonho112.github.io/svyder/reference/glance.svyder.md),
[`print.svyder()`](https://joonho112.github.io/svyder/reference/print.svyder.md),
[`summary.svyder()`](https://joonho112.github.io/svyder/reference/summary.svyder.md),
[`tidy.svyder()`](https://joonho112.github.io/svyder/reference/tidy.svyder.md)

## Examples

``` r
is.svyder(list())
#> [1] FALSE

data(nsece_demo)
result <- der_compute(
  nsece_demo$draws,
  y = nsece_demo$y, X = nsece_demo$X,
  group = nsece_demo$group, weights = nsece_demo$weights,
  psu = nsece_demo$psu, family = "binomial",
  sigma_theta = nsece_demo$sigma_theta,
  param_types = nsece_demo$param_types
)
is.svyder(result)
#> [1] TRUE
```
