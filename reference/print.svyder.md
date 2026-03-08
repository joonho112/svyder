# Print a svyder Object

Displays a concise summary of DER diagnostic results, including the DER
range, classification threshold, number of flagged parameters, and
correction status.

## Usage

``` r
# S3 method for class 'svyder'
print(x, n = 10, digits = 3, ...)
```

## Arguments

- x:

  A `svyder` object.

- n:

  Maximum number of flagged parameters to display (default 10).

- digits:

  Number of decimal places for DER values (default 3).

- ...:

  Ignored.

## Value

Invisibly returns `x`.

## See also

[`summary.svyder()`](https://joonho112.github.io/svyder/reference/summary.svyder.md)
for detailed classification output,
[`tidy.svyder()`](https://joonho112.github.io/svyder/reference/tidy.svyder.md)
for a tidy data frame.

Other svyder-methods:
[`as.matrix.svyder()`](https://joonho112.github.io/svyder/reference/as.matrix.svyder.md),
[`glance.svyder()`](https://joonho112.github.io/svyder/reference/glance.svyder.md),
[`is.svyder()`](https://joonho112.github.io/svyder/reference/is.svyder.md),
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
print(result)
#> svyder diagnostic (54 parameters)
#>   Family: binomial | N = 6785 | J = 51
#>   DER range: [0.235, 5.315]
#>   Threshold (tau): 1.20
#>   Flagged: 30 / 54 (55.6%)
#> 
#>   Flagged parameters:
#>     beta[2]              DER = 2.687  [I-a] -> CORRECT
#>     theta[1]             DER = 3.384  [II] -> CORRECT
#>     theta[4]             DER = 2.212  [II] -> CORRECT
#>     theta[5]             DER = 1.571  [II] -> CORRECT
#>     theta[6]             DER = 2.103  [II] -> CORRECT
#>     theta[7]             DER = 2.241  [II] -> CORRECT
#>     theta[9]             DER = 5.315  [II] -> CORRECT
#>     theta[11]            DER = 2.653  [II] -> CORRECT
#>     theta[15]            DER = 1.573  [II] -> CORRECT
#>     theta[18]            DER = 4.022  [II] -> CORRECT
#>     ... and 20 more
#> 
#>   Correction applied: 30 parameter(s) rescaled
#>   Compute time: 0.100 sec
```
