# Build the clustered score outer product matrix

Constructs J_cluster = sum_g s_g s_g^T where s_g is the score vector
summed within PSU g.

## Usage

``` r
.build_J_cluster(X, r, psu, group, p, J)
```

## Arguments

- X:

  Design matrix (N x p).

- r:

  Weighted residuals (length N).

- psu:

  Integer PSU indicator (1:G), length N.

- group:

  Integer group indicator (1:J), length N.

- p:

  Number of fixed-effect parameters.

- J:

  Number of groups.

## Value

Symmetric (p+J) x (p+J) matrix.
