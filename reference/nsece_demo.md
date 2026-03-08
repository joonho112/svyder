# Synthetic NSECE-Like Survey Data

A synthetic dataset mimicking the NSECE 2019 survey structure for
demonstration of the DER diagnostic pipeline. Contains N = 6785
observations across J = 51 states with unequal survey weights, clustered
PSU structure, and three fixed-effect covariates (intercept,
within-cluster poverty, between-cluster tiered reimbursement policy).

## Usage

``` r
nsece_demo
```

## Format

A list with components:

- draws:

  Matrix of posterior draws (4000 x 54), columns 1:3 are fixed effects
  (beta), columns 4:54 are random effects (theta).

- y:

  Binary outcome vector (length 6785).

- X:

  Design matrix (6785 x 3) with columns: intercept, poverty_cwc
  (group-mean centered), tiered_reim (binary policy).

- group:

  Integer state group indicator (1 to 51).

- weights:

  Survey weights (positive, length 6785). Log-normal distributed,
  normalized within state.

- psu:

  PSU indicators (integer, length 6785).

- param_types:

  Character vector of length 3:
  `c("fe_between", "fe_within", "fe_between")`.

- family:

  Model family: `"binomial"`.

- sigma_theta:

  Random effect SD (0.66).

- N:

  Number of observations (6785).

- J:

  Number of groups (51).

- p:

  Number of fixed effects (3).

## Source

Synthetic data generated to mimic NSECE 2019 structure. See
`data-raw/generate_nsece_demo.R`.

## Examples

``` r
data(nsece_demo)
str(nsece_demo, max.level = 1)
#> List of 12
#>  $ draws      : num [1:4000, 1:54] 0.3321 0.0921 0.1227 -0.0035 0.2898 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ y          : int [1:6785] 0 0 0 0 0 1 0 1 0 0 ...
#>  $ X          : num [1:6785, 1:3] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ group      : int [1:6785] 1 1 1 1 1 1 1 1 1 1 ...
#>  $ weights    : num [1:6785] 0.327 6.848 0.269 0.272 0.937 ...
#>  $ psu        : num [1:6785] 1 1 1 1 1 1 1 1 1 1 ...
#>  $ param_types: chr [1:3] "fe_between" "fe_within" "fe_between"
#>  $ family     : chr "binomial"
#>  $ sigma_theta: num 0.66
#>  $ N          : num 6785
#>  $ J          : int 51
#>  $ p          : int 3
```
