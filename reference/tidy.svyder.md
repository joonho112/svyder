# Tidy a svyder Object

Returns a one-row-per-parameter data frame with DER diagnostics in a
tidy format compatible with broom conventions. Includes posterior
summaries, classification tier, and correction scale factor.

## Usage

``` r
tidy.svyder(x, ...)
```

## Arguments

- x:

  A `svyder` object.

- ...:

  Additional arguments (unused).

## Value

A `data.frame` with one row per parameter and columns:

- term:

  Parameter name.

- estimate:

  Posterior mean (from original draws).

- std.error:

  Posterior standard deviation (from original draws).

- der:

  Design Effect Ratio.

- tier:

  Three-tier classification (if classified).

- action:

  Action label: `"CORRECT"` or `"retain"` (if classified).

- flagged:

  Logical; whether the parameter is flagged for correction (if
  classified).

- scale_factor:

  Cholesky scale factor applied to this parameter.

## See also

[`glance.svyder()`](https://joonho112.github.io/svyder/reference/glance.svyder.md)
for model-level summaries,
[`print.svyder()`](https://joonho112.github.io/svyder/reference/print.svyder.md)
for console output.

Other svyder-methods:
[`as.matrix.svyder()`](https://joonho112.github.io/svyder/reference/as.matrix.svyder.md),
[`glance.svyder()`](https://joonho112.github.io/svyder/reference/glance.svyder.md),
[`is.svyder()`](https://joonho112.github.io/svyder/reference/is.svyder.md),
[`print.svyder()`](https://joonho112.github.io/svyder/reference/print.svyder.md),
[`summary.svyder()`](https://joonho112.github.io/svyder/reference/summary.svyder.md)

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
head(tidy.svyder(result))
#>              term   estimate  std.error       der tier  action flagged
#> beta[1]   beta[1]  0.2498445 0.14735061 0.2617696  I-b  retain   FALSE
#> beta[2]   beta[2] -0.1494408 0.02604573 2.6868825  I-a CORRECT    TRUE
#> beta[3]   beta[3]  0.1610078 0.20239837 0.3427294  I-b  retain   FALSE
#> theta[1] theta[1] -0.2281087 0.40501395 3.3838218   II CORRECT    TRUE
#> theta[2] theta[2]  0.7500209 0.42644395 0.6758803   II  retain   FALSE
#> theta[3] theta[3]  1.0856688 0.43304778 1.1187639   II  retain   FALSE
#>          scale_factor
#> beta[1]      1.000000
#> beta[2]      1.639171
#> beta[3]      1.000000
#> theta[1]     1.839517
#> theta[2]     1.000000
#> theta[3]     1.000000
```
