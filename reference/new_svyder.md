# Low-Level Constructor for svyder Objects

Creates a svyder object from pre-computed DER pipeline results. This is
the internal constructor; users should call
[`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md)
or
[`der_diagnose()`](https://joonho112.github.io/svyder/reference/der_diagnose.md)
instead.

## Usage

``` r
new_svyder(
  der,
  params,
  H_obs,
  J_cluster,
  V_sand,
  sigma_mcmc,
  deff_j,
  B_j,
  classification,
  tau,
  corrected_draws,
  scale_factors,
  original_draws,
  call,
  family,
  n_obs,
  n_groups,
  compute_time
)
```

## Arguments

- der:

  Named numeric vector of Design Effect Ratios.

- params:

  Character vector of parameter names.

- H_obs:

  Numeric matrix, observed information matrix (d x d).

- J_cluster:

  Numeric matrix, clustered score outer product (d x d).

- V_sand:

  Numeric matrix, sandwich variance (d x d).

- sigma_mcmc:

  Numeric matrix, MCMC posterior covariance (d x d).

- deff_j:

  Numeric vector of cluster-level design effects (length J).

- B_j:

  Numeric vector of shrinkage factors (length J).

- classification:

  Data frame with per-parameter classification results.

- tau:

  Numeric scalar, flagging threshold for DER.

- corrected_draws:

  Numeric matrix, corrected posterior draws (M x d).

- scale_factors:

  Numeric vector of Cholesky scale factors (length d).

- original_draws:

  Numeric matrix, original posterior draws (M x d).

- call:

  The matched call that produced this object.

- family:

  Character string, model family (e.g., "binomial", "gaussian").

- n_obs:

  Integer, number of observations.

- n_groups:

  Integer, number of groups/clusters.

- compute_time:

  Numeric scalar, computation time in seconds.

## Value

A list of class `"svyder"`.
