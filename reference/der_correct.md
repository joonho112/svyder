# Apply Selective Correction to Flagged Parameters

Rescales the posterior draws of parameters flagged by
[`der_classify()`](https://joonho112.github.io/svyder/reference/der_classify.md)
so that their covariance matches the sandwich variance target, leaving
unflagged parameters untouched. The correction preserves posterior
means.

## Usage

``` r
der_correct(x, method = c("block_cholesky", "marginal"))
```

## Arguments

- x:

  A `svyder` object with classification (from
  [`der_classify()`](https://joonho112.github.io/svyder/reference/der_classify.md)).

- method:

  `"block_cholesky"` (default) or `"marginal"`.

## Value

A `svyder` object with `corrected_draws`, `scale_factors` (marginal SD
ratios of flagged parameters), and `correction_method` populated.

## Details

With `method = "block_cholesky"` (the paper's Algorithm 1), the flagged
block F is transformed jointly: writing \\\Sigma_F = L_1 L_1^\top\\ and
\\V_F = L_2 L_2^\top\\ for the MCMC and sandwich covariance of the
flagged block, each draw is mapped through \\\phi^\* = \hat\phi + L_2
L_1^{-1} (\phi - \hat\phi)\\, so the corrected block has covariance
exactly \\V_F\\. With a single flagged parameter this reduces to scaling
by \\\sqrt{V\_{ii}/\Sigma\_{ii}}\\.

With `method = "marginal"`, each flagged parameter is rescaled
independently by \\\sqrt{V\_{ii}/\Sigma\_{ii}}\\. This matches marginal
variances only and is provided for transparency and for reproducing
results computed this way; it is not a Cholesky correction when more
than one parameter is flagged.

## References

Lee, J., Williams, M. R., & Savitsky, T. D. (2026). Design Effect Ratios
for Bayesian Survey Models: A Diagnostic Framework for Identifying
Survey-Sensitive Parameters. *Journal of Survey Statistics and
Methodology*. Submitted.

Savitsky, T. D., & Williams, M. R. (2022). Pseudo Bayesian mixed models
under informative sampling. *Journal of Official Statistics*, 38(4),
901–928.
[doi:10.2478/jos-2022-0039](https://doi.org/10.2478/jos-2022-0039)

Williams, M. R., & Savitsky, T. D. (2021). Uncertainty estimation for
pseudo-Bayesian inference under complex sampling. *International
Statistical Review*, 89(1), 72–107.
[doi:10.1111/insr.12376](https://doi.org/10.1111/insr.12376)

Higham, N. J. (2002). Computing the nearest correlation matrix—a problem
from finance. *IMA Journal of Numerical Analysis*, 22(3), 329–343.
[doi:10.1093/imanum/22.3.329](https://doi.org/10.1093/imanum/22.3.329)

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
  cluster = nsece_demo$psu, family = "binomial",
  sigma_theta = nsece_demo$sigma_theta,
  param_types = nsece_demo$param_types
)
result <- der_classify(result, tau = 1.2, verbose = FALSE)
result <- der_correct(result)
```
