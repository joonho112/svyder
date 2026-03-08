# Extract Survey Design Information

Extracts weights, cluster identifiers, and strata from survey design
objects created by the survey package. The extracted information can be
passed directly to
[`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md).

## Usage

``` r
extract_design(design, ...)

# S3 method for class 'survey.design2'
extract_design(design, ...)

# Default S3 method
extract_design(design, ...)
```

## Arguments

- design:

  A survey design object (class `survey.design2`).

- ...:

  Additional arguments passed to methods.

## Value

A list with components:

- weights:

  Numeric vector of survey weights.

- cluster:

  Factor or integer vector of PSU/cluster identifiers.

- strata:

  Data frame of strata variables (may be `NULL` for unstratified
  designs).

## See also

[`extract_draws()`](https://joonho112.github.io/svyder/reference/extract_draws.md)
for extracting posterior draws,
[`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md)
for computing DER.

Other extraction:
[`extract_draws.brmsfit()`](https://joonho112.github.io/svyder/reference/extract_draws.md)
