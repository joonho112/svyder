# Verify Theoretical DER Predictions

Compares empirical DER values against theoretical predictions from the
decomposition theorems. This function checks how well the closed-form
approximations match the numerically computed DER.

## Usage

``` r
der_theorem_check(x)
```

## Arguments

- x:

  A `svyder` object.

## Value

A `data.frame` with columns: `param`, `param_type`, `der_empirical`,
`der_theorem1` (for FE), `der_theorem2` (for RE), `relative_error`,
`theorem_used`. If the conservation law is applicable, the result has a
`"conservation_law"` attribute.

## Details

For fixed effects (Theorem 1):

- `fe_within`: DER \\\approx\\ DEFF.

- `fe_between`: DER \\\approx\\ DEFF \\\cdot\\ (1 - B).

For random effects (Theorem 2):

- DER \\\approx\\ B \\\cdot\\ DEFF \\\cdot\\ kappa(J).

Also reports the conservation-law diagnostic (Corollary 1) when
applicable: DER_mu + DER_theta_cond \\\approx\\ DEFF (balanced
intercept-only case). The conditional random-effect term is reported
separately from the empirical marginal mean of the fitted random-effect
DERs.

## References

Lee, J., Williams, M. R., & Savitsky, T. D. (2026). Design Effect Ratios
for Bayesian Survey Models: A Diagnostic Framework for Identifying
Survey-Sensitive Parameters. *Journal of Survey Statistics and
Methodology*. Submitted.

## See also

[`der_decompose()`](https://joonho112.github.io/svyder/reference/der_decompose.md)
for the full decomposition.

Other analysis:
[`der_compare()`](https://joonho112.github.io/svyder/reference/der_compare.md),
[`der_decompose()`](https://joonho112.github.io/svyder/reference/der_decompose.md),
[`der_sensitivity()`](https://joonho112.github.io/svyder/reference/der_sensitivity.md)

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
thm <- der_theorem_check(result)
head(thm)
#>      param param_type der_empirical der_theorem1 der_theorem2 relative_error
#> 1  beta[1] fe_between     0.2626885    0.3781171           NA     0.43941210
#> 2  beta[2]  fe_within     2.6894416    2.5952698           NA     0.03501540
#> 3  beta[3] fe_between     0.3430527    0.3781171           NA     0.10221260
#> 4 theta[1]         re     3.3858856           NA     2.122360     0.37317436
#> 5 theta[2]         re     0.6759045           NA     1.396910     1.06672753
#> 6 theta[3]         re     1.1182861           NA     1.048017     0.06283619
#>          theorem_used
#> 1 Theorem 1 (between)
#> 2  Theorem 1 (within)
#> 3 Theorem 1 (between)
#> 4      Theorem 2 (RE)
#> 5      Theorem 2 (RE)
#> 6      Theorem 2 (RE)
```
