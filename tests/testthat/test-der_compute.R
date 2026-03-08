# test-der_compute.R
# Tests for the core der_compute() function


# ============================================================================
# Helper: run der_compute on a fixture list
# ============================================================================
.run_der_compute <- function(fix) {
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)
  der_compute(
    draws_all,
    y            = fix$y,
    X            = fix$X,
    group        = fix$group,
    weights      = fix$w,
    psu          = fix$psu,
    family       = fix$family,
    sigma_theta  = fix$sigma_theta_hat,
    sigma_e      = fix$sigma_e,
    beta_prior_sd = fix$beta_prior_sd,
    param_types  = fix$param_types

  )
}


# ============================================================================
# Basic functionality
# ============================================================================

test_that("der_compute returns svyder class", {
  fix <- make_balanced_gaussian()
  result <- .run_der_compute(fix)
  expect_s3_class(result, "svyder")
})

test_that("DER values are positive", {
  fix <- make_balanced_gaussian()
  result <- .run_der_compute(fix)
  expect_true(all(result$der > 0))
})

test_that("DER length matches parameters", {
  fix <- make_balanced_gaussian()
  result <- .run_der_compute(fix)
  expected_d <- fix$p + fix$J
  expect_equal(length(result$der), expected_d)
  expect_equal(length(result$params), expected_d)
})

test_that("pipe-friendly: returns svyder", {
  fix <- make_balanced_gaussian()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)
  result <- der_compute(
    draws_all,
    y = fix$y, X = fix$X, group = fix$group,
    weights = fix$w, psu = fix$psu,
    family = "gaussian", sigma_theta = fix$sigma_theta_hat,
    sigma_e = fix$sigma_e, beta_prior_sd = fix$beta_prior_sd,
    param_types = fix$param_types
  )
  expect_s3_class(result, "svyder")
  # Can pipe into der_classify
  result2 <- der_classify(result, tau = 1.2, verbose = FALSE)
  expect_s3_class(result2, "svyder")
})

test_that("errors on unknown class", {
  expect_error(
    der_compute("not_a_matrix"),
    "does not know how to handle class"
  )
  expect_error(
    der_compute(data.frame(x = 1:10)),
    "does not know how to handle class"
  )
})


# ============================================================================
# Regression test: matches standalone compute_der()
# ============================================================================

test_that("matches standalone compute_der()", {
  standalone_path <- file.path(
    "/Users/joonholee/Documents/2026-03-03_Design Effect Ratios Paper",
    "dev/scripts/standalone/compute_der.R"
  )
  skip_if_not(file.exists(standalone_path),
              "Standalone compute_der.R not found")

  source(standalone_path, local = TRUE)

  # Use binomial fixture (standalone only supports binomial)
  fix <- make_unbalanced_binomial()

  # --- Run standalone ---
  standalone_result <- compute_der(
    beta_hat         = fix$beta_hat,
    theta_hat        = fix$theta_hat,
    sigma_theta_hat  = fix$sigma_theta_hat,
    draws_beta       = fix$draws_beta,
    draws_theta      = fix$draws_theta,
    y                = fix$y,
    X                = fix$X,
    group            = fix$group,
    w                = fix$w,
    psu              = fix$psu,
    beta_prior_sd    = fix$beta_prior_sd,
    tau              = 1.2,
    param_types      = fix$param_types
  )

  # --- Run svyder ---
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)
  svyder_result <- der_compute(
    draws_all,
    y            = fix$y,
    X            = fix$X,
    group        = fix$group,
    weights      = fix$w,
    psu          = fix$psu,
    family       = "binomial",
    sigma_theta  = fix$sigma_theta_hat,
    beta_prior_sd = fix$beta_prior_sd,
    param_types  = fix$param_types
  )

  # DER values should be identical (same math path)
  # Note: der_compute uses colMeans(draws) as point estimates while
  # standalone uses fixed beta_hat/theta_hat. We need to account for this.
  # The standalone uses the true values as point_est, while der_compute
  # uses colMeans of draws. Since draws are generated from the true values,
  # the colMeans will be close but not identical. So we compare with a
  # small tolerance rather than exact equality.
  expect_equal(
    as.numeric(svyder_result$der),
    as.numeric(standalone_result$der),
    tolerance = 0.05,
    info = "DER values should closely match standalone"
  )

  # Sandwich matrices should also be close
  expect_equal(
    diag(svyder_result$V_sand),
    diag(standalone_result$V_sand),
    tolerance = 0.05
  )
})


# ============================================================================
# Gaussian fixture
# ============================================================================

test_that("balanced gaussian works end-to-end", {
  fix <- make_balanced_gaussian()
  result <- .run_der_compute(fix)

  # Should have correct dimensions
  expect_equal(result$n_obs, fix$N)
  expect_equal(result$n_groups, fix$J)
  expect_equal(result$family, "gaussian")

  # DER for equal weights should be positive and finite
  # Some random effects may have small DER due to strong shrinkage
  expect_true(all(result$der > 0))
  expect_true(all(is.finite(result$der)))
})


# ============================================================================
# Binomial fixture
# ============================================================================

test_that("binomial fixture works end-to-end", {
  fix <- make_unbalanced_binomial()
  result <- .run_der_compute(fix)

  expect_s3_class(result, "svyder")
  expect_equal(result$family, "binomial")
  expect_equal(result$n_obs, fix$N)
  expect_equal(result$n_groups, fix$J)
  expect_equal(length(result$der), fix$d)
  expect_true(all(result$der > 0))
})


# ============================================================================
# Minimal J=2 fixture
# ============================================================================

test_that("minimal J=2 fixture works", {
  fix <- make_minimal_j2()
  result <- .run_der_compute(fix)

  expect_s3_class(result, "svyder")
  expect_equal(result$n_groups, 2L)
  expect_equal(length(result$deff_j), 2L)
  expect_equal(length(result$B_j), 2L)
})


# ============================================================================
# Stored components
# ============================================================================

test_that("all matrix components have correct dimensions", {
  fix <- make_unbalanced_binomial()
  result <- .run_der_compute(fix)
  d <- fix$d

  expect_equal(dim(result$H_obs), c(d, d))
  expect_equal(dim(result$J_cluster), c(d, d))
  expect_equal(dim(result$V_sand), c(d, d))
  expect_equal(dim(result$sigma_mcmc), c(d, d))
  expect_equal(dim(result$original_draws), c(fix$M, d))
})

test_that("sigma_mcmc is symmetric positive definite", {
  fix <- make_balanced_gaussian()
  result <- .run_der_compute(fix)

  # Symmetric
  expect_equal(result$sigma_mcmc, t(result$sigma_mcmc))

  # Positive definite (all eigenvalues positive)
  eigenvals <- eigen(result$sigma_mcmc, only.values = TRUE)$values
  expect_true(all(eigenvals > 0))
})

test_that("V_sand is symmetric", {
  fix <- make_balanced_gaussian()
  result <- .run_der_compute(fix)
  expect_equal(result$V_sand, t(result$V_sand))
})
