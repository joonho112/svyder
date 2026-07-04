# Declared cluster universe (cluster_strata): selected-but-empty clusters
# must enter the meat as zero score-total clusters rather than silently
# vanishing from the design target.

test_that(".build_J_cluster matches a hand-rolled meat with an empty cluster", {
  set.seed(42)
  # 2 strata (= groups), 3 selected clusters each; cluster 6 has NO rows.
  N <- 10L
  p <- 2L
  J <- 2L
  X <- cbind(1, rnorm(N))
  r <- rnorm(N)
  cluster <- c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L)
  group   <- c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L)
  strata  <- group
  cluster_strata <- c(1L, 1L, 1L, 2L, 2L, 2L)   # cluster 6: stratum 2, empty

  J_built <- svyder:::.build_J_cluster(
    X, r, cluster, group, p, J,
    strata = strata, cluster_strata = cluster_strata,
    center = TRUE, df_correct = TRUE
  )

  # Hand-rolled: universe of G = 6 clusters, zero row for cluster 6.
  G <- 6L
  d <- p + J
  T_mat <- matrix(0, G, d)
  for (i in seq_len(N)) {
    g <- cluster[i]
    T_mat[g, 1:p] <- T_mat[g, 1:p] + r[i] * X[i, ]
    T_mat[g, p + group[i]] <- T_mat[g, p + group[i]] + r[i]
  }
  J_man <- matrix(0, d, d)
  for (h in 1:2) {
    idx <- which(cluster_strata == h)
    C_h <- length(idx)
    T_h <- sweep(T_mat[idx, , drop = FALSE], 2,
                 colMeans(T_mat[idx, , drop = FALSE]))
    J_man <- J_man + C_h / (C_h - 1) * crossprod(T_h)
  }
  J_man <- (J_man + t(J_man)) / 2

  expect_equal(J_built, J_man, tolerance = 1e-12)

  # And it must DIFFER from the dropped-cluster meat (the pre-fix behavior):
  J_dropped <- svyder:::.build_J_cluster(
    X, r, cluster, group, p, J,
    strata = strata,
    center = TRUE, df_correct = TRUE
  )
  expect_gt(max(abs(J_built - J_dropped)), 1e-8)
})

test_that(".build_J_cluster validates cluster_strata inputs", {
  X <- cbind(1, rnorm(4)); r <- rnorm(4)
  cluster <- c(1L, 1L, 2L, 2L); group <- c(1L, 1L, 2L, 2L)

  expect_error(
    svyder:::.build_J_cluster(X, r, cluster, group, 2L, 2L,
                              strata = NULL,
                              cluster_strata = c(1L, 2L, 2L)),
    "requires 'strata'"
  )
  expect_error(
    svyder:::.build_J_cluster(X, r, cluster, group, 2L, 2L,
                              strata = c(1L, 1L, 1L, 1L),
                              cluster_strata = c(1L, 2L, 2L)),
    "disagrees"
  )
})

test_that("der_compute records and uses a declared cluster universe", {
  set.seed(7)
  J <- 3L; n_per <- 30L; N <- J * n_per; p <- 2L
  group  <- rep(seq_len(J), each = n_per)
  X      <- cbind(intercept = 1, x = rnorm(N))
  theta  <- c(-0.4, 0.1, 0.5)
  y      <- rbinom(N, 1L, plogis(X %*% c(0.2, 0.5) + theta[group]))
  w      <- runif(N, 0.5, 2)
  draws  <- matrix(rnorm(200 * (p + J), sd = 0.3), 200, p + J)

  # 3 selected PSUs per group-stratum; PSU 9 (group 3) is empty.
  psu_of_unit <- rep(rep(1:2, length.out = n_per), J) +
    3L * (group - 1L)                       # units land only in PSUs 1,2 mod 3
  cluster_strata <- rep(seq_len(J), each = 3L)

  res <- der_compute(
    draws, y = y, X = X, group = group, weights = w,
    cluster = psu_of_unit, strata = group,
    cluster_strata = cluster_strata,
    family = "binomial", sigma_theta = 0.5,
    param_types = c("fe_within", "fe_within")
  )
  expect_s3_class(res, "svyder")
  expect_identical(res$target$cluster_n, 9L)
  expect_identical(res$target$cluster_n_empty, 3L)   # PSUs 3, 6, 9 empty
  expect_true(res$target$cluster_universe_declared)

  # Without the universe: fewer clusters, different meat.
  res0 <- der_compute(
    draws, y = y, X = X, group = group, weights = w,
    cluster = psu_of_unit, strata = group,
    family = "binomial", sigma_theta = 0.5,
    param_types = c("fe_within", "fe_within")
  )
  expect_identical(res0$target$cluster_n, 6L)
  expect_gt(max(abs(res$J_cluster - res0$J_cluster)), 1e-10)
})
