# Classify Parameters by Design Sensitivity

Assigns each parameter to a design-sensitivity tier and flags those
whose DER exceeds the threshold `tau`. The three-tier classification:

- **Tier I-a** (`fe_within`): Survey-dominated parameters.

- **Tier I-b** (`fe_between`): Protected between-cluster parameters.

- **Tier II** (`re`): Protected random effects.

## Usage

``` r
der_classify(x, tau = 1.2, verbose = TRUE)
```

## Arguments

- x:

  A `svyder` object from
  [`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md).

- tau:

  Threshold (default 1.2). Parameters with DER \> tau are flagged.

- verbose:

  Print classification summary (default `TRUE`).

## Value

A `svyder` object with updated `classification` and `tau` fields.

## Details

Parameters with DER \> `tau` (strict inequality) are flagged for
correction, regardless of tier.

## See also

[`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md)
for computing DER values,
[`der_correct()`](https://joonho112.github.io/svyder/reference/der_correct.md)
for applying corrections,
[`der_diagnose()`](https://joonho112.github.io/svyder/reference/der_diagnose.md)
for the all-in-one pipeline.

Other core-pipeline:
[`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md),
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
result <- der_classify(result, tau = 1.2)
#> DER Classification (tau = 1.20)
#>   Total parameters: 54
#>   Flagged: 30 (55.6%)
#>   Flagged parameters:
#>     beta[2]: DER = 2.687 [I-a] -> CORRECT
#>     theta[1]: DER = 3.384 [II] -> CORRECT
#>     theta[4]: DER = 2.212 [II] -> CORRECT
#>     theta[5]: DER = 1.571 [II] -> CORRECT
#>     theta[6]: DER = 2.103 [II] -> CORRECT
#>     theta[7]: DER = 2.241 [II] -> CORRECT
#>     theta[9]: DER = 5.315 [II] -> CORRECT
#>     theta[11]: DER = 2.653 [II] -> CORRECT
#>     theta[15]: DER = 1.573 [II] -> CORRECT
#>     theta[18]: DER = 4.022 [II] -> CORRECT
#>     theta[19]: DER = 1.992 [II] -> CORRECT
#>     theta[20]: DER = 2.477 [II] -> CORRECT
#>     theta[21]: DER = 1.790 [II] -> CORRECT
#>     theta[23]: DER = 1.311 [II] -> CORRECT
#>     theta[27]: DER = 2.326 [II] -> CORRECT
#>     theta[29]: DER = 1.632 [II] -> CORRECT
#>     theta[30]: DER = 2.972 [II] -> CORRECT
#>     theta[31]: DER = 1.290 [II] -> CORRECT
#>     theta[32]: DER = 1.713 [II] -> CORRECT
#>     theta[34]: DER = 2.842 [II] -> CORRECT
#>     theta[35]: DER = 1.259 [II] -> CORRECT
#>     theta[36]: DER = 2.277 [II] -> CORRECT
#>     theta[39]: DER = 2.222 [II] -> CORRECT
#>     theta[40]: DER = 1.567 [II] -> CORRECT
#>     theta[41]: DER = 1.818 [II] -> CORRECT
#>     theta[43]: DER = 2.300 [II] -> CORRECT
#>     theta[44]: DER = 3.233 [II] -> CORRECT
#>     theta[45]: DER = 1.515 [II] -> CORRECT
#>     theta[46]: DER = 1.228 [II] -> CORRECT
#>     theta[47]: DER = 2.141 [II] -> CORRECT
```
