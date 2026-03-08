# Apply Selective Cholesky Correction

For each parameter flagged by
[`der_classify()`](https://joonho112.github.io/svyder/reference/der_classify.md),
rescales the posterior draws so that marginal variance matches the
sandwich variance estimate. The correction preserves the posterior mean:
draws are centered, scaled by `sqrt(V_sand[i,i] / sigma_mcmc[i,i])`,
then re-centered.

## Usage

``` r
der_correct(x, method = "cholesky")
```

## Arguments

- x:

  A `svyder` object with classification (from
  [`der_classify()`](https://joonho112.github.io/svyder/reference/der_classify.md)).

- method:

  Correction method (default `"cholesky"`). Currently only `"cholesky"`
  is supported.

## Value

A `svyder` object with `corrected_draws`, `scale_factors`, and
`original_draws` populated.

## Details

Unflagged parameters retain their original draws without any
modification.

## See also

[`der_classify()`](https://joonho112.github.io/svyder/reference/der_classify.md)
for flagging parameters,
[`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md)
for computing DER,
[`as.matrix.svyder()`](https://joonho112.github.io/svyder/reference/as.matrix.svyder.md)
for extracting corrected draws.

Other core-pipeline:
[`der_classify()`](https://joonho112.github.io/svyder/reference/der_classify.md),
[`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md),
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
result <- der_classify(result, tau = 1.2, verbose = FALSE)
result <- der_correct(result)
```
