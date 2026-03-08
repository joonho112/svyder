# Map rstanarm parameter names to svyder convention

Converts rstanarm-style parameter names to the svyder naming convention:

- `(Intercept)` -\> `beta[1]`

- `x1` -\> `beta[2]`

- `b[(Intercept) group:level1]` -\> `theta[1]`

- `b[(Intercept) group:level2]` -\> `theta[2]`

## Usage

``` r
.map_rstanarm_params(param_names)
```

## Arguments

- param_names:

  Character vector of rstanarm parameter names.

## Value

Data frame with columns: original_name, svyder_name, type.

## Details

Excludes variance components (`Sigma[...]`, `sigma`) and log-posterior.
