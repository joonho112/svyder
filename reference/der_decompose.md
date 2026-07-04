# Decompose DER into Components

Decomposes each parameter's DER into its constituent factors. For fixed
effects, the protection factor \\R_k\\ is computed structurally from the
working-weighted within/between split of covariate \\k\\'s identifying
variation: \$\$R_k = \frac{\sum_j B_j\\ a_j \bar x\_{jk}^2} {\sum_j
(W\_{jkk} + a_j \bar x\_{jk}^2)},\$\$ where \\a_j\\ is the working
information of group \\j\\, \\\bar x\_{jk}\\ the working-weighted group
mean of covariate \\k\\, and \\W\_{jkk}\\ its within-group working sum
of squares. A pure within-group covariate has \\R_k = 0\\; a pure
between-group covariate has \\R_k\\ equal to the information-weighted
average shrinkage factor. The prediction is \\DER_k \approx DEFF \cdot
(1 - R_k)\\ with the global Kish DEFF.

## Usage

``` r
der_decompose(x)
```

## Arguments

- x:

  A `svyder` object from svyder \>= 0.2.0 (the stored `$data` slots are
  required).

## Value

A `data.frame` with one row per parameter and columns:

- param:

  Parameter name.

- param_type:

  `"fe_within"`, `"fe_between"`, or `"re"`.

- der:

  Observed Design Effect Ratio.

- deff_used:

  Design effect entering the prediction: the global Kish DEFF for fixed
  effects, the per-group \\\mathrm{DEFF}\_j\\ for random effects.

- B_used:

  Shrinkage factor \\B_j\\ (random-effect rows only; `NA` for fixed
  effects).

- R_k:

  Structural between-group share of covariate \\k\\'s identifying
  variation (fixed-effect rows only; `NA` for random effects).

- kappa:

  Finite-group correction \\\kappa_j\\ (random-effect rows only; `NA`
  for fixed effects).

- der_predicted:

  Theorem-based prediction: \\\mathrm{DEFF}\\ (1 - R_k)\\ for fixed
  effects (Theorem 1), \\B_j\\\mathrm{DEFF}\_j\\\kappa_j\\ for random
  effects (Theorem 2).

## Details

For random effects the per-group prediction \\DER_j \approx B_j \cdot
DEFF_j \cdot \kappa_j\\ is used (matching
[`der_theorem_check()`](https://joonho112.github.io/svyder/reference/der_theorem_check.md)),
not a pooled mean.

Because `der_predicted` is structural, comparing it with the observed
`der` is a genuine check of the decomposition theorems' mechanism,
within their stated model class.

The prediction is row-specific: `der_predicted` uses Theorem 1
(\\\mathrm{DEFF}\cdot(1 - R_k)\\) for fixed-effect rows and Theorem 2
(\\B_j\cdot\mathrm{DEFF}\_j\cdot\kappa_j\\) for random-effect rows.
Consequently the component columns are populated per row type: `B_used`
and `kappa` are `NA` for fixed effects (whose prediction uses only DEFF
and \\R_k\\), and `R_k` is `NA` for random effects (whose prediction
uses only \\B_j\\, \\\mathrm{DEFF}\_j\\ and \\\kappa_j\\).

## References

Lee, J., Williams, M. R., & Savitsky, T. D. (2026). Design Effect Ratios
for Bayesian Survey Models: A Diagnostic Framework for Identifying
Survey-Sensitive Parameters. *Journal of Survey Statistics and
Methodology*. Submitted.

Kish, L. (1965). *Survey Sampling*. John Wiley & Sons.

Fay, R. E., & Herriot, R. A. (1979). Estimates of income for small
places: An application of James–Stein procedures to census data.
*Journal of the American Statistical Association*, 74(366), 269–277.
[doi:10.1080/01621459.1979.10482505](https://doi.org/10.1080/01621459.1979.10482505)

## See also

[`der_theorem_check()`](https://joonho112.github.io/svyder/reference/der_theorem_check.md)
for theorem-level comparison,
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
  cluster = nsece_demo$psu, family = "binomial",
  sigma_theta = nsece_demo$sigma_theta,
  param_types = nsece_demo$param_types
)
decomp <- der_decompose(result)
head(decomp)
#>      param param_type       der deff_used    B_used        R_k     kappa
#> 1  beta[1] fe_between 0.2626885  2.727828        NA 0.93273657        NA
#> 2  beta[2]  fe_within 2.6894416  2.727828        NA 0.01430246        NA
#> 3  beta[3] fe_between 0.3430527  2.727828        NA 0.92282567        NA
#> 4 theta[1]         re 3.3858856  3.479991 0.6441524         NA 0.9467869
#> 5 theta[2]         re 0.6759045  2.417686 0.6072076         NA 0.9515495
#> 6 theta[3]         re 1.1182861  1.878577 0.5847477         NA 0.9540496
#>   der_predicted
#> 1     0.1834831
#> 2     2.6888138
#> 3     0.2105183
#> 4     2.1223599
#> 5     1.3969105
#> 6     1.0480173
```
