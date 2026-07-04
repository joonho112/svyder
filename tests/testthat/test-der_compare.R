# test-der_compare.R
# Tests for DER comparison across clustering definitions

test_that("compare returns correct columns", {
  fix <- make_balanced_gaussian()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)

  clusters <- list(
    state = fix$psu,
    same  = fix$psu
  )

  cmp <- der_compare(
    draws_all,
    clusters = clusters,
    y = fix$y,
    X = fix$X,
    group = fix$group,
    weights = fix$w,
    family = fix$family,
    sigma_theta = fix$sigma_theta_hat,
    sigma_e = fix$sigma_e,
    beta_prior_sd = fix$beta_prior_sd,
    param_types = fix$param_types
  )

  expect_s3_class(cmp, "data.frame")
  expected_cols <- c("param", "cluster_name", "der")
  expect_true(all(expected_cols %in% names(cmp)))
  expect_equal(nrow(cmp), 2 * fix$d)  # 2 clusterings * d params
})


test_that("different clusterings give different DER", {
  fix <- make_unbalanced_binomial()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)

  # Create a different clustering: split each group into 2 PSUs
  psu_fine <- numeric(fix$N)
  psu_counter <- 0L
  for (j in seq_len(fix$J)) {
    idx_j <- which(fix$group == j)
    half  <- ceiling(length(idx_j) / 2)
    psu_counter <- psu_counter + 1L
    psu_fine[idx_j[seq_len(half)]] <- psu_counter
    psu_counter <- psu_counter + 1L
    psu_fine[idx_j[seq(half + 1, length(idx_j))]] <- psu_counter
  }

  clusters <- list(
    coarse = fix$psu,
    fine   = as.integer(psu_fine)
  )

  cmp <- der_compare(
    draws_all,
    clusters = clusters,
    y = fix$y,
    X = fix$X,
    group = fix$group,
    weights = fix$w,
    family = fix$family,
    sigma_theta = fix$sigma_theta_hat,
    beta_prior_sd = fix$beta_prior_sd,
    param_types = fix$param_types
  )

  # DER values should differ between the two clusterings
  coarse_der <- cmp$der[cmp$cluster_name == "coarse"]
  fine_der   <- cmp$der[cmp$cluster_name == "fine"]
  expect_false(all(abs(coarse_der - fine_der) < 1e-10))
})


test_that("same clustering gives identical DER", {
  fix <- make_balanced_gaussian()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)

  clusters <- list(
    run1 = fix$psu,
    run2 = fix$psu
  )

  cmp <- der_compare(
    draws_all,
    clusters = clusters,
    y = fix$y,
    X = fix$X,
    group = fix$group,
    weights = fix$w,
    family = fix$family,
    sigma_theta = fix$sigma_theta_hat,
    sigma_e = fix$sigma_e,
    beta_prior_sd = fix$beta_prior_sd,
    param_types = fix$param_types
  )

  run1_der <- cmp$der[cmp$cluster_name == "run1"]
  run2_der <- cmp$der[cmp$cluster_name == "run2"]
  expect_equal(run1_der, run2_der)
})


test_that("compare requires named clusters list", {
  fix <- make_balanced_gaussian()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)

  # Unnamed list should error
  expect_error(
    der_compare(
      draws_all,
      clusters = list(fix$psu),
      y = fix$y,
      X = fix$X,
      group = fix$group,
      weights = fix$w,
      family = fix$family,
      sigma_theta = fix$sigma_theta_hat,
      sigma_e = fix$sigma_e,
      beta_prior_sd = fix$beta_prior_sd,
      param_types = fix$param_types
    ),
    "named list"
  )
})


test_that("compare works with single clustering", {
  fix <- make_balanced_gaussian()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)

  clusters <- list(default = fix$psu)

  cmp <- der_compare(
    draws_all,
    clusters = clusters,
    y = fix$y,
    X = fix$X,
    group = fix$group,
    weights = fix$w,
    family = fix$family,
    sigma_theta = fix$sigma_theta_hat,
    sigma_e = fix$sigma_e,
    beta_prior_sd = fix$beta_prior_sd,
    param_types = fix$param_types
  )

  expect_equal(nrow(cmp), fix$d)
  expect_true(all(cmp$cluster_name == "default"))
})


test_that("svyder fast path equals matrix-path full recomputation", {
  fix <- make_unbalanced_binomial()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)

  # A finer clustering: split each group into 2 clusters
  psu_fine <- integer(fix$N)
  counter <- 0L
  for (j in seq_len(fix$J)) {
    idx_j <- which(fix$group == j)
    half  <- ceiling(length(idx_j) / 2)
    counter <- counter + 1L
    psu_fine[idx_j[seq_len(half)]] <- counter
    counter <- counter + 1L
    psu_fine[idx_j[seq(half + 1, length(idx_j))]] <- counter
  }

  clusters <- list(
    model_group = fix$group,
    fine_psu    = psu_fine
  )

  # Fast path: reuse the stored data slots of a svyder object
  sv <- der_compute(
    draws_all,
    y = fix$y, X = fix$X, group = fix$group,
    weights = fix$w, cluster = fix$psu,
    family = "binomial", sigma_theta = fix$sigma_theta_hat,
    beta_prior_sd = fix$beta_prior_sd, param_types = fix$param_types
  )
  cmp_fast <- der_compare(sv, clusters = clusters)

  # Matrix path: full recomputation per target
  cmp_full <- der_compare(
    draws_all,
    clusters = clusters,
    y = fix$y, X = fix$X, group = fix$group,
    weights = fix$w,
    family = "binomial", sigma_theta = fix$sigma_theta_hat,
    beta_prior_sd = fix$beta_prior_sd, param_types = fix$param_types
  )

  expect_equal(cmp_fast$param, cmp_full$param)
  expect_equal(cmp_fast$cluster_name, cmp_full$cluster_name)
  expect_equal(cmp_fast$der, cmp_full$der, tolerance = 1e-12)

  # The model-group target from the fast path also matches a direct
  # der_compute() with cluster = group
  sv_group <- der_compute(
    draws_all,
    y = fix$y, X = fix$X, group = fix$group,
    weights = fix$w, cluster = fix$group,
    family = "binomial", sigma_theta = fix$sigma_theta_hat,
    beta_prior_sd = fix$beta_prior_sd, param_types = fix$param_types
  )
  fast_group <- cmp_fast$der[cmp_fast$cluster_name == "model_group"]
  expect_equal(fast_group, as.numeric(sv_group$der), tolerance = 1e-12)
})


test_that("svyder fast path preserves declared cluster universe for stored target", {
  set.seed(7)
  J <- 3L; n_per <- 30L; N <- J * n_per; p <- 2L
  group <- rep(seq_len(J), each = n_per)
  X <- cbind(intercept = 1, x = rnorm(N))
  theta <- c(-0.4, 0.1, 0.5)
  y <- rbinom(N, 1L, plogis(X %*% c(0.2, 0.5) + theta[group]))
  w <- runif(N, 0.5, 2)
  draws <- matrix(rnorm(200 * (p + J), sd = 0.3), 200, p + J)

  # 3 selected PSUs per group-stratum; PSUs 3, 6, and 9 are empty.
  psu_of_unit <- rep(rep(1:2, length.out = n_per), J) +
    3L * (group - 1L)
  cluster_strata <- rep(seq_len(J), each = 3L)

  res <- der_compute(
    draws, y = y, X = X, group = group, weights = w,
    cluster = psu_of_unit, strata = group,
    cluster_strata = cluster_strata,
    family = "binomial", sigma_theta = 0.5,
    param_types = c("fe_within", "fe_within")
  )

  cmp_same <- der_compare(res, clusters = list(same = psu_of_unit))
  expect_equal(cmp_same$der, as.numeric(res$der), tolerance = 1e-12)

  cmp_group <- suppressWarnings(
    der_compare(res, clusters = list(model_group = group))
  )
  direct_group <- suppressWarnings(
    der_compute(
      draws, y = y, X = X, group = group, weights = w,
      cluster = group, strata = group,
      family = "binomial", sigma_theta = 0.5,
      param_types = c("fe_within", "fe_within")
    )
  )
  expect_equal(cmp_group$der, as.numeric(direct_group$der), tolerance = 1e-12)
})
