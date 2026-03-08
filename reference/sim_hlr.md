# Simulated Hierarchical Linear Regression Data

A small balanced Gaussian hierarchical model dataset for quick testing
and demonstration. Contains J = 10 groups with n_j = 20 observations
each (N = 200 total), equal weights (DEFF = 1), and two fixed-effect
covariates (intercept + within-cluster covariate).

## Usage

``` r
sim_hlr
```

## Format

A list with components:

- draws:

  Matrix of posterior draws (4000 x 12), columns 1:2 are fixed effects
  (beta), columns 3:12 are random effects (theta).

- y:

  Continuous outcome vector (length 200).

- X:

  Design matrix (200 x 2) with columns: intercept, x_within.

- group:

  Integer group indicator (1 to 10).

- weights:

  Survey weights (all 1.0, length 200).

- psu:

  PSU indicators (same as group).

- param_types:

  Character vector of length 2: `c("fe_between", "fe_within")`.

- family:

  Model family: `"gaussian"`.

- sigma_theta:

  Random effect SD (0.5).

- sigma_e:

  Residual SD (1.0).

- N:

  Number of observations (200).

- J:

  Number of groups (10).

- p:

  Number of fixed effects (2).

- B_ref:

  Analytical shrinkage factor (5/6).

- deff_ref:

  Reference design effect (1.0).

## Source

Synthetic data. See `data-raw/generate_sim_hlr.R`.

## Details

With equal weights the design effect is 1.0, so DER values should be
close to 1.0 across all parameters. This dataset is useful for verifying
that the pipeline correctly identifies the absence of design effects.

## Examples

``` r
data(sim_hlr)
str(sim_hlr, max.level = 1)
#> List of 15
#>  $ draws      : num [1:4000, 1:12] 2.09 2.05 2.11 1.98 1.93 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ y          : num [1:200] 1.173 1.73 2.698 1.094 0.951 ...
#>  $ X          : num [1:200, 1:2] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ group      : int [1:200] 1 1 1 1 1 1 1 1 1 1 ...
#>  $ weights    : num [1:200] 1 1 1 1 1 1 1 1 1 1 ...
#>  $ psu        : int [1:200] 1 1 1 1 1 1 1 1 1 1 ...
#>  $ param_types: chr [1:2] "fe_between" "fe_within"
#>  $ family     : chr "gaussian"
#>  $ sigma_theta: num 0.5
#>  $ sigma_e    : num 1
#>  $ N          : int 200
#>  $ J          : int 10
#>  $ p          : int 2
#>  $ B_ref      : num 0.833
#>  $ deff_ref   : num 1
```
