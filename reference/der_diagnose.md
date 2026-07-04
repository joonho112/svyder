# All-in-One DER Diagnostic

Runs the full DER pipeline in a single call:
[`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md)
to obtain DER values,
[`der_classify()`](https://joonho112.github.io/svyder/reference/der_classify.md)
to assign tiers and flag parameters, and optionally
[`der_correct()`](https://joonho112.github.io/svyder/reference/der_correct.md)
to apply selective correction. This is the recommended entry point for
most users.

## Usage

``` r
der_diagnose(
  x,
  ...,
  tau = 1.2,
  correct = TRUE,
  method = c("block_cholesky", "marginal")
)
```

## Arguments

- x:

  A draws matrix (S x d), brmsfit, CmdStanMCMC, or stanreg object. For
  the matrix method, columns 1:p are fixed effects and columns
  (p+1):(p+J) are random effects.

- ...:

  Additional arguments passed to methods.

- tau:

  Classification threshold (default 1.2).

- correct:

  Apply correction to flagged parameters (default `TRUE`).

- method:

  Correction method passed to
  [`der_correct()`](https://joonho112.github.io/svyder/reference/der_correct.md)
  (default `"block_cholesky"`).

## Value

A fully processed `svyder` object with DER values, classification, and
(optionally) corrected draws.

## References

Lee, J., Williams, M. R., & Savitsky, T. D. (2026). Design Effect Ratios
for Bayesian Survey Models: A Diagnostic Framework for Identifying
Survey-Sensitive Parameters. *Journal of Survey Statistics and
Methodology*. Submitted.

## See also

[`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md),
[`der_classify()`](https://joonho112.github.io/svyder/reference/der_classify.md),
[`der_correct()`](https://joonho112.github.io/svyder/reference/der_correct.md)
for the individual pipeline steps.
[`tidy.svyder()`](https://joonho112.github.io/svyder/reference/tidy.svyder.md)
and
[`glance.svyder()`](https://joonho112.github.io/svyder/reference/glance.svyder.md)
for tidy summaries.

Other core-pipeline:
[`der_classify()`](https://joonho112.github.io/svyder/reference/der_classify.md),
[`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md),
[`der_correct()`](https://joonho112.github.io/svyder/reference/der_correct.md)

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
print(result)
#> svyder diagnostic (54 parameters)
#>   Family: binomial | N = 6785 | J = 51
#>   Target: 644 cluster(s) | meat: centered, DF-corrected | weights: unit_mean
#>   DER range: [0.235, 5.308]
#>   Tier III (DER undefined): sigma_theta
#>   Threshold (tau): 1.20
#>   Flagged: 30 / 54 (55.6%)
#> 
#>   Flagged parameters:
#>     beta[2]              DER = 2.689  [I-a] -> CORRECT
#>     theta[1]             DER = 3.386  [II] -> CORRECT
#>     theta[4]             DER = 2.210  [II] -> CORRECT
#>     theta[5]             DER = 1.568  [II] -> CORRECT
#>     theta[6]             DER = 2.104  [II] -> CORRECT
#>     theta[7]             DER = 2.232  [II] -> CORRECT
#>     theta[9]             DER = 5.308  [II] -> CORRECT
#>     theta[11]            DER = 2.655  [II] -> CORRECT
#>     theta[15]            DER = 1.576  [II] -> CORRECT
#>     theta[18]            DER = 4.025  [II] -> CORRECT
#>     ... and 20 more
#> 
#>   Correction applied: 30 parameter(s) rescaled
#>   Compute time: 0.057 sec
summary(result)
#>  param_name param_type       der tier                 tier_label flagged
#>     beta[1] fe_between 0.2626885  I-b        Protected (between)   FALSE
#>     beta[2]  fe_within 2.6894416  I-a           Survey-dominated    TRUE
#>     beta[3] fe_between 0.3430527  I-b        Protected (between)   FALSE
#>    theta[1]         re 3.3858856   II Protected (random effects)    TRUE
#>    theta[2]         re 0.6759045   II Protected (random effects)   FALSE
#>    theta[3]         re 1.1182861   II Protected (random effects)   FALSE
#>    theta[4]         re 2.2098098   II Protected (random effects)    TRUE
#>    theta[5]         re 1.5677648   II Protected (random effects)    TRUE
#>    theta[6]         re 2.1043742   II Protected (random effects)    TRUE
#>    theta[7]         re 2.2322332   II Protected (random effects)    TRUE
#>    theta[8]         re 0.5997084   II Protected (random effects)   FALSE
#>    theta[9]         re 5.3077337   II Protected (random effects)    TRUE
#>   theta[10]         re 0.3239584   II Protected (random effects)   FALSE
#>   theta[11]         re 2.6550681   II Protected (random effects)    TRUE
#>   theta[12]         re 0.5000038   II Protected (random effects)   FALSE
#>   theta[13]         re 0.2354213   II Protected (random effects)   FALSE
#>   theta[14]         re 0.5529524   II Protected (random effects)   FALSE
#>   theta[15]         re 1.5756844   II Protected (random effects)    TRUE
#>   theta[16]         re 0.6803656   II Protected (random effects)   FALSE
#>   theta[17]         re 0.8177403   II Protected (random effects)   FALSE
#>   theta[18]         re 4.0246235   II Protected (random effects)    TRUE
#>   theta[19]         re 1.9922108   II Protected (random effects)    TRUE
#>   theta[20]         re 2.4691218   II Protected (random effects)    TRUE
#>   theta[21]         re 1.7903889   II Protected (random effects)    TRUE
#>   theta[22]         re 0.8754189   II Protected (random effects)   FALSE
#>   theta[23]         re 1.3105172   II Protected (random effects)    TRUE
#>   theta[24]         re 0.9748128   II Protected (random effects)   FALSE
#>   theta[25]         re 0.7496725   II Protected (random effects)   FALSE
#>   theta[26]         re 1.1578875   II Protected (random effects)   FALSE
#>   theta[27]         re 2.3241799   II Protected (random effects)    TRUE
#>   theta[28]         re 1.1860962   II Protected (random effects)   FALSE
#>   theta[29]         re 1.6306621   II Protected (random effects)    TRUE
#>   theta[30]         re 2.9678609   II Protected (random effects)    TRUE
#>   theta[31]         re 1.2876481   II Protected (random effects)    TRUE
#>   theta[32]         re 1.7154958   II Protected (random effects)    TRUE
#>   theta[33]         re 0.7052294   II Protected (random effects)   FALSE
#>   theta[34]         re 2.8436195   II Protected (random effects)    TRUE
#>   theta[35]         re 1.2581949   II Protected (random effects)    TRUE
#>   theta[36]         re 2.2803658   II Protected (random effects)    TRUE
#>   theta[37]         re 0.7436107   II Protected (random effects)   FALSE
#>   theta[38]         re 1.1182233   II Protected (random effects)   FALSE
#>   theta[39]         re 2.2223976   II Protected (random effects)    TRUE
#>   theta[40]         re 1.5668657   II Protected (random effects)    TRUE
#>   theta[41]         re 1.8209951   II Protected (random effects)    TRUE
#>   theta[42]         re 0.8845594   II Protected (random effects)   FALSE
#>   theta[43]         re 2.3013513   II Protected (random effects)    TRUE
#>   theta[44]         re 3.2167380   II Protected (random effects)    TRUE
#>   theta[45]         re 1.5048423   II Protected (random effects)    TRUE
#>   theta[46]         re 1.2253161   II Protected (random effects)    TRUE
#>   theta[47]         re 2.1341967   II Protected (random effects)    TRUE
#>   theta[48]         re 0.8883745   II Protected (random effects)   FALSE
#>   theta[49]         re 0.7717870   II Protected (random effects)   FALSE
#>   theta[50]         re 0.8962989   II Protected (random effects)   FALSE
#>   theta[51]         re 0.7652428   II Protected (random effects)   FALSE
#>   action
#>   retain
#>  CORRECT
#>   retain
#>  CORRECT
#>   retain
#>   retain
#>  CORRECT
#>  CORRECT
#>  CORRECT
#>  CORRECT
#>   retain
#>  CORRECT
#>   retain
#>  CORRECT
#>   retain
#>   retain
#>   retain
#>  CORRECT
#>   retain
#>   retain
#>  CORRECT
#>  CORRECT
#>  CORRECT
#>  CORRECT
#>   retain
#>  CORRECT
#>   retain
#>   retain
#>   retain
#>  CORRECT
#>   retain
#>  CORRECT
#>  CORRECT
#>  CORRECT
#>  CORRECT
#>   retain
#>  CORRECT
#>  CORRECT
#>  CORRECT
#>   retain
#>   retain
#>  CORRECT
#>  CORRECT
#>  CORRECT
#>   retain
#>  CORRECT
#>  CORRECT
#>  CORRECT
#>  CORRECT
#>  CORRECT
#>   retain
#>   retain
#>   retain
#>   retain
```
