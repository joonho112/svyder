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
