# Decompose DER into Components

Decomposes each parameter's DER into its constituent factors: Kish DEFF,
shrinkage factor B, protection factor R_k, and finite-J correction
kappa. This decomposition reveals why each parameter has its observed
DER value.

## Usage

``` r
der_decompose(x)
```

## Arguments

- x:

  A `svyder` object.

## Value

A `data.frame` with columns: `param`, `param_type`, `der`, `deff_mean`,
`B_mean`, `R_k`, `kappa`, `der_predicted`.

## Details

For fixed effects:

- `fe_within`: DER \\\approx\\ DEFF \\\cdot\\ (1 - R_k), where R_k
  \\\approx\\ 0.

- `fe_between`: DER \\\approx\\ DEFF \\\cdot\\ (1 - R_k), where R_k
  \\\approx\\ B.

For random effects:

- DER \\\approx\\ B \\\cdot\\ DEFF \\\cdot\\ kappa(J).

## See also

[`der_theorem_check()`](https://joonho112.github.io/svyder/reference/der_theorem_check.md)
for verifying theoretical predictions,
[`plot.svyder()`](https://joonho112.github.io/svyder/reference/plot.svyder.md)
with `type = "decomposition"` for visualization.

Other analysis:
[`der_compare()`](https://joonho112.github.io/svyder/reference/der_compare.md),
[`der_sensitivity()`](https://joonho112.github.io/svyder/reference/der_sensitivity.md),
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
decomp <- der_decompose(result)
head(decomp)
#>      param param_type       der deff_mean    B_mean       R_k     kappa
#> 1  beta[1] fe_between 0.2617696   2.59527 0.8543053 0.8991359        NA
#> 2  beta[2]  fe_within 2.6868825   2.59527 0.8543053 0.0000000        NA
#> 3  beta[3] fe_between 0.3427294   2.59527 0.8543053 0.8679407        NA
#> 4 theta[1]         re 3.3838218   2.59527 0.8543053        NA 0.8213246
#> 5 theta[2]         re 0.6758803   2.59527 0.8543053        NA 0.8213246
#> 6 theta[3]         re 1.1187639   2.59527 0.8543053        NA 0.8213246
#>   der_predicted
#> 1     0.2617696
#> 2     2.5952698
#> 3     0.3427294
#> 4     1.8210021
#> 5     1.8210021
#> 6     1.8210021
```
