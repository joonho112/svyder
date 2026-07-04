# Build the clustered score outer product matrix (meat)

Constructs the meat matrix from cluster-level score totals. The general
(default) form is the stratified, centered, DF-corrected estimator \$\$J
= \sum_h \frac{C_h}{C_h - 1} \sum_c (t\_{hc} - \bar t_h)(t\_{hc} - \bar
t_h)^\top,\$\$ where \\t\_{hc}\\ is the score total of cluster \\c\\ in
stratum \\h\\, \\C_h\\ the number of clusters in stratum \\h\\, and
\\\bar t_h\\ the stratum mean of cluster totals. Centering matters at
the posterior-mean plug-in point because the likelihood score totals do
not sum to zero there (the prior gradient does not vanish).

## Usage

``` r
.build_J_cluster(
  X,
  r,
  cluster,
  group,
  p,
  J,
  strata = NULL,
  cluster_strata = NULL,
  center = TRUE,
  df_correct = TRUE
)
```

## Arguments

- X:

  Design matrix (N x p).

- r:

  Weighted score residuals (length N).

- cluster:

  Integer cluster indicator (1:G), length N. This is the sandwich
  aggregation unit (design PSU or model group).

- group:

  Integer group indicator (1:J), length N.

- p:

  Number of fixed-effect parameters.

- J:

  Number of groups.

- strata:

  Optional integer stratum indicator (1:H), length N. `NULL` treats the
  design as a single stratum.

- cluster_strata:

  Optional integer vector giving the stratum of every cluster in the
  declared design universe, indexed by cluster id (length G, the
  universe size). Supplying it declares the cluster universe explicitly:
  cluster ids in `1:G` that appear in no observation – e.g., selected
  PSUs whose second-stage Poisson sample is empty – remain in the meat
  as zero score-total clusters, entering the stratum centering and the
  C_h/(C_h-1) correction. Without it, the universe is inferred from the
  observed ids, which drops empty clusters. Requires `strata`.

- center:

  Center cluster score totals within strata (default TRUE).

- df_correct:

  Apply the C_h / (C_h - 1) stratum correction (default TRUE).

## Value

Symmetric (p+J) x (p+J) matrix.

## Details

With `strata = NULL`, `center = FALSE`, `df_correct = FALSE` this
reduces to the uncentered single-stratum form \\J = \sum_c t_c
t_c^\top\\ used by the v1 pipeline (kept for exact reproduction of
previously published results).

Strata containing a single cluster cannot be centered; their uncentered
contribution is used with a warning.
