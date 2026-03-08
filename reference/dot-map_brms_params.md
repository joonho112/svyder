# Map brms parameter names to svyder convention

Converts brms-style parameter names to the svyder naming convention:

- `b_Intercept` -\> `beta[1]`

- `b_x1` -\> `beta[2]`

- `r_group[level1,Intercept]` -\> `theta[1]`

- `r_group[level2,Intercept]` -\> `theta[2]`

## Usage

``` r
.map_brms_params(param_names)
```

## Arguments

- param_names:

  Character vector of brms parameter names.

## Value

Data frame with columns:

- original_name:

  The original brms parameter name.

- svyder_name:

  The mapped svyder parameter name.

- type:

  Either "fixed" or "random".

## Details

Excludes variance components (`sd_`, `cor_`), residual SD (`sigma`), and
log-posterior (`lp__`, `lprior`).
