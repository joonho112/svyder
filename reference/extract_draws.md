# Extract Posterior Draws from Model Objects

S3 generic for extracting posterior draws from fitted model objects.
Each method returns a standardized list with a draws matrix and optional
parameter metadata, which can then be passed to
[`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md).

## Usage

``` r
# S3 method for class 'brmsfit'
extract_draws(x, ..., pars = NULL)

# S3 method for class 'CmdStanMCMC'
extract_draws(x, ..., pars = NULL)

# S3 method for class 'draws_matrix'
extract_draws(x, ...)

# S3 method for class 'draws_df'
extract_draws(x, ...)

# S3 method for class 'draws_array'
extract_draws(x, ...)

# S3 method for class 'draws_list'
extract_draws(x, ...)

# S3 method for class 'draws_rvars'
extract_draws(x, ...)

# S3 method for class 'stanreg'
extract_draws(x, ..., pars = NULL)

extract_draws(x, ...)

# S3 method for class 'matrix'
extract_draws(x, ...)

# Default S3 method
extract_draws(x, ...)
```

## Arguments

- x:

  A fitted model object (matrix, brmsfit, CmdStanMCMC, stanreg,
  draws_matrix, or draws_df).

- ...:

  Additional arguments passed to methods.

- pars:

  Character vector of parameter names to extract from the stanreg
  object. If `NULL` (default), extracts all fixed effect and group-level
  parameters, excluding variance components and log-posterior.

## Value

A list with components:

- draws:

  Numeric matrix of posterior draws (S x d), where S is the number of
  posterior samples and d is the number of parameters.

- param_info:

  (Optional) Data frame with parameter metadata including original names
  and mapped names.

## Details

Methods are provided for `matrix` (identity), `brmsfit`, `CmdStanMCMC`,
`stanreg`, `draws_matrix`, `draws_df`, `draws_array`, `draws_list`, and
`draws_rvars`.

## References

Lee, J., Williams, M. R., & Savitsky, T. D. (2026). Design Effect Ratios
for Bayesian Survey Models: A Diagnostic Framework for Identifying
Survey-Sensitive Parameters. *Journal of Survey Statistics and
Methodology*. Submitted.

Rabe-Hesketh, S., & Skrondal, A. (2006). Multilevel modelling of complex
survey data. *Journal of the Royal Statistical Society: Series A*,
169(4), 805–827.
[doi:10.1111/j.1467-985X.2006.00426.x](https://doi.org/10.1111/j.1467-985X.2006.00426.x)

## See also

[`extract_design()`](https://joonho112.github.io/svyder/reference/extract_design.md)
for extracting survey design information,
[`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md)
for computing DER from draws.

Other backends:
[`extract_design()`](https://joonho112.github.io/svyder/reference/extract_design.md)

## Examples

``` r
# Matrix method is the identity
m <- matrix(rnorm(100), nrow = 20)
result <- extract_draws(m)
identical(result$draws, m)
#> [1] TRUE
```
