# Safely invert a matrix with nearPD fallback

Attempts `solve(H_obs)`. If that fails (singular or near-singular),
falls back to
[`Matrix::nearPD`](https://rdrr.io/pkg/Matrix/man/nearPD.html) to find
the nearest positive-definite matrix, then inverts that.

## Usage

``` r
.safe_invert(H_obs)
```

## Arguments

- H_obs:

  Square numeric matrix.

## Value

The inverse of H_obs (or its nearest PD approximation).
