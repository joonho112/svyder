# Build the observed information matrix

Constructs the (p+J) x (p+J) block matrix H_obs from working weights.

## Usage

``` r
.build_H_obs(X, v, group, J, p, beta_prior_prec, theta_prior_prec)
```

## Arguments

- X:

  Design matrix (N x p).

- v:

  Working weights (length N).

- group:

  Integer group indicator (1:J), length N.

- J:

  Number of groups.

- p:

  Number of fixed-effect parameters.

- beta_prior_prec:

  Scalar prior precision for beta.

- theta_prior_prec:

  Scalar prior precision for theta.

## Value

Symmetric (p+J) x (p+J) matrix.
