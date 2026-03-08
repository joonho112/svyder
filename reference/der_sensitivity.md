# Sensitivity Analysis Across Threshold Values

Evaluates how the number of flagged parameters changes across a range of
threshold values `tau`. Useful for assessing the robustness of
classification results to the choice of threshold.

## Usage

``` r
der_sensitivity(x, tau_range = seq(0.8, 2, by = 0.1))
```

## Arguments

- x:

  A `svyder` object.

- tau_range:

  Numeric vector of threshold values to evaluate. Default:
  `seq(0.8, 2.0, by = 0.1)`.

## Value

A `data.frame` with columns: `tau`, `n_flagged`, `pct_flagged`, and
`flagged_params` (a list-column of character vectors naming the flagged
parameters at each threshold).

## See also

[`der_classify()`](https://joonho112.github.io/svyder/reference/der_classify.md)
for classification at a single threshold.

Other analysis:
[`der_compare()`](https://joonho112.github.io/svyder/reference/der_compare.md),
[`der_decompose()`](https://joonho112.github.io/svyder/reference/der_decompose.md),
[`der_theorem_check()`](https://joonho112.github.io/svyder/reference/der_theorem_check.md)

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
sens <- der_sensitivity(result)
sens[, c("tau", "n_flagged", "pct_flagged")]
#>    tau n_flagged pct_flagged
#> 1  0.8        40   0.7407407
#> 2  0.9        35   0.6481481
#> 3  1.0        34   0.6296296
#> 4  1.1        34   0.6296296
#> 5  1.2        30   0.5555556
#> 6  1.3        27   0.5000000
#> 7  1.4        26   0.4814815
#> 8  1.5        26   0.4814815
#> 9  1.6        22   0.4074074
#> 10 1.7        21   0.3888889
#> 11 1.8        19   0.3518519
#> 12 1.9        18   0.3333333
#> 13 2.0        17   0.3148148
```
