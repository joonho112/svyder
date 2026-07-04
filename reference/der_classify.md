# Classify Parameters by Design Sensitivity

Assigns each estimable parameter to a design-sensitivity tier and flags
those whose DER exceeds the threshold `tau`. The estimable-parameter
classification:

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

## Tier DER signatures

Each tier has a closed-form DER signature under the decomposition
theorems (with \\B\\ the shrinkage factor, \\J\\ the number of groups):

- **Tier I-a** (within-group fixed effects): \\\mathrm{DER} =
  \mathrm{DEFF}\\ – fully exposed to the design and the primary
  correction candidates.

- **Tier I-b** (between-group fixed effects): \\\mathrm{DER} =
  \mathrm{DEFF}\\(1 - B)\\ – shielded by shrinkage and usually below the
  threshold.

- **Tier II** (random effects): \\\mathrm{DER} =
  B\cdot\mathrm{DEFF}\cdot\kappa(J)\\, with the finite-group correction
  \\\kappa(J) = 1 - 1/\[J(1-B)+B\]\\.

- **Tier III** (hyperparameters, e.g. \\\sigma\_\theta\\): excluded by
  construction because they are outside the declared location-block DER
  target. Reported in `$excluded`, never flagged.

## References

Lee, J., Williams, M. R., & Savitsky, T. D. (2026). Design Effect Ratios
for Bayesian Survey Models: A Diagnostic Framework for Identifying
Survey-Sensitive Parameters. *Journal of Survey Statistics and
Methodology*. Submitted.

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
  cluster = nsece_demo$psu, family = "binomial",
  sigma_theta = nsece_demo$sigma_theta,
  param_types = nsece_demo$param_types
)
result <- der_classify(result, tau = 1.2)
#> DER Classification (tau = 1.20)
#>   Total parameters: 54
#>   Tier III excluded (DER undefined): sigma_theta
#>   Flagged: 30 (55.6%)
#>   Flagged parameters:
#>     beta[2]: DER = 2.689 [I-a] -> CORRECT
#>     theta[1]: DER = 3.386 [II] -> CORRECT
#>     theta[4]: DER = 2.210 [II] -> CORRECT
#>     theta[5]: DER = 1.568 [II] -> CORRECT
#>     theta[6]: DER = 2.104 [II] -> CORRECT
#>     theta[7]: DER = 2.232 [II] -> CORRECT
#>     theta[9]: DER = 5.308 [II] -> CORRECT
#>     theta[11]: DER = 2.655 [II] -> CORRECT
#>     theta[15]: DER = 1.576 [II] -> CORRECT
#>     theta[18]: DER = 4.025 [II] -> CORRECT
#>     theta[19]: DER = 1.992 [II] -> CORRECT
#>     theta[20]: DER = 2.469 [II] -> CORRECT
#>     theta[21]: DER = 1.790 [II] -> CORRECT
#>     theta[23]: DER = 1.311 [II] -> CORRECT
#>     theta[27]: DER = 2.324 [II] -> CORRECT
#>     theta[29]: DER = 1.631 [II] -> CORRECT
#>     theta[30]: DER = 2.968 [II] -> CORRECT
#>     theta[31]: DER = 1.288 [II] -> CORRECT
#>     theta[32]: DER = 1.715 [II] -> CORRECT
#>     theta[34]: DER = 2.844 [II] -> CORRECT
#>     theta[35]: DER = 1.258 [II] -> CORRECT
#>     theta[36]: DER = 2.280 [II] -> CORRECT
#>     theta[39]: DER = 2.222 [II] -> CORRECT
#>     theta[40]: DER = 1.567 [II] -> CORRECT
#>     theta[41]: DER = 1.821 [II] -> CORRECT
#>     theta[43]: DER = 2.301 [II] -> CORRECT
#>     theta[44]: DER = 3.217 [II] -> CORRECT
#>     theta[45]: DER = 1.505 [II] -> CORRECT
#>     theta[46]: DER = 1.225 [II] -> CORRECT
#>     theta[47]: DER = 2.134 [II] -> CORRECT
```
