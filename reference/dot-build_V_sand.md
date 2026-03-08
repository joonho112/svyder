# Build the sandwich variance estimator

Computes V_sand = H_obs_inv %*% J_cluster %*% H_obs_inv and symmetrizes.

## Usage

``` r
.build_V_sand(H_obs_inv, J_cluster)
```

## Arguments

- H_obs_inv:

  Inverse of the observed information matrix (d x d).

- J_cluster:

  Clustered score outer product matrix (d x d).

## Value

Symmetric (d x d) sandwich variance matrix.
