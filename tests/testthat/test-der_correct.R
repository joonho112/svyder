# test-der_correct.R
# Tests for der_correct()


# ============================================================================
# Helper: full pipeline through correction
# ============================================================================
.full_pipeline <- function(fix, tau = 1.2) {
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)
  result <- der_compute(
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
  result <- der_classify(result, tau = tau, verbose = FALSE)
  der_correct(result)
}


# ============================================================================
# Scale factor tests
# ============================================================================

test_that("scale factors = sqrt(DER) for flagged", {
  fix <- make_unbalanced_binomial()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)
  result <- der_compute(
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

  # Force some parameters to be flagged
  result$der[1] <- 2.0
  result$der[2] <- 1.5
  result <- der_classify(result, tau = 1.2, verbose = FALSE)
  result <- der_correct(result)

  flagged_idx <- which(result$classification$flagged)
  for (i in flagged_idx) {
    expected_sf <- unname(sqrt(diag(result$V_sand)[i] / diag(result$sigma_mcmc)[i]))
    expect_equal(result$scale_factors[i], expected_sf, tolerance = 1e-10)
  }
})

test_that("scale factors = 1.0 for unflagged", {
  fix <- make_balanced_gaussian()
  result <- .full_pipeline(fix, tau = 1.2)

  unflagged_idx <- which(!result$classification$flagged)
  for (i in unflagged_idx) {
    expect_equal(result$scale_factors[i], 1.0)
  }
})


# ============================================================================
# Mean preservation
# ============================================================================

test_that("mean preservation", {
  fix <- make_unbalanced_binomial()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)
  result <- der_compute(
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

  # Force flagging
  result$der[1] <- 2.5
  result <- der_classify(result, tau = 1.2, verbose = FALSE)
  result <- der_correct(result)

  # Column means of corrected draws should match original draws
  original_means  <- colMeans(result$original_draws)
  corrected_means <- colMeans(result$corrected_draws)

  expect_equal(corrected_means, original_means, tolerance = 1e-10)
})


# ============================================================================
# Unflagged draws bitwise identical
# ============================================================================

test_that("unflagged draws bitwise identical", {
  fix <- make_unbalanced_binomial()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)
  result <- der_compute(
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

  # Force one parameter flagged, rest unflagged
  result$der[1] <- 2.0
  result <- der_classify(result, tau = 1.2, verbose = FALSE)
  result <- der_correct(result)

  unflagged_idx <- which(!result$classification$flagged)
  for (i in unflagged_idx) {
    # Bitwise identical: use identical() not expect_equal()
    expect_identical(
      result$corrected_draws[, i],
      result$original_draws[, i],
      info = paste("Column", i, "should be bitwise identical")
    )
  }
})


# ============================================================================
# Corrected variance matches target
# ============================================================================

test_that("corrected variance matches target for flagged params", {
  fix <- make_unbalanced_binomial()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)
  result <- der_compute(
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

  # Force flagging
  result$der[1] <- 2.0
  result <- der_classify(result, tau = 1.2, verbose = FALSE)
  result <- der_correct(result)

  flagged_idx <- which(result$classification$flagged)
  for (i in flagged_idx) {
    corrected_var <- var(result$corrected_draws[, i])
    target_var    <- diag(result$V_sand)[i]
    # The corrected variance should be close to the sandwich variance
    # (not exactly equal because var() uses n-1 denominator)
    expect_equal(corrected_var, target_var, tolerance = 0.05)
  }
})


# ============================================================================
# No flagged -> draws unchanged
# ============================================================================

test_that("no flagged -> draws unchanged", {
  fix <- make_balanced_gaussian()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)
  result <- der_compute(
    draws_all,
    y            = fix$y,
    X            = fix$X,
    group        = fix$group,
    weights      = fix$w,
    psu          = fix$psu,
    family       = "gaussian",
    sigma_theta  = fix$sigma_theta_hat,
    sigma_e      = fix$sigma_e,
    beta_prior_sd = fix$beta_prior_sd,
    param_types  = fix$param_types
  )

  # Set all DER below threshold
  result$der[] <- 0.9
  result <- der_classify(result, tau = 1.2, verbose = FALSE)
  result <- der_correct(result)

  expect_identical(result$corrected_draws, result$original_draws)
  expect_true(all(result$scale_factors == 1.0))
})


# ============================================================================
# Full pipeline
# ============================================================================

test_that("full pipeline: der_compute |> der_classify |> der_correct", {
  fix <- make_unbalanced_binomial()
  result <- .full_pipeline(fix, tau = 1.2)

  expect_s3_class(result, "svyder")
  expect_true(!is.null(result$corrected_draws))
  expect_true(!is.null(result$original_draws))
  expect_equal(ncol(result$corrected_draws), fix$d)
  expect_equal(nrow(result$corrected_draws), fix$M)
  expect_equal(length(result$scale_factors), fix$d)
})


# ============================================================================
# as.matrix.svyder
# ============================================================================

test_that("as.matrix.svyder extracts corrected draws", {
  fix <- make_unbalanced_binomial()
  result <- .full_pipeline(fix, tau = 1.2)

  mat <- as.matrix(result)
  expect_true(is.matrix(mat))
  expect_equal(dim(mat), dim(result$corrected_draws))
  expect_identical(mat, result$corrected_draws)
})

test_that("as.matrix.svyder falls back to original draws when no correction", {
  fix <- make_balanced_gaussian()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)
  result <- der_compute(
    draws_all,
    y            = fix$y,
    X            = fix$X,
    group        = fix$group,
    weights      = fix$w,
    psu          = fix$psu,
    family       = "gaussian",
    sigma_theta  = fix$sigma_theta_hat,
    sigma_e      = fix$sigma_e,
    beta_prior_sd = fix$beta_prior_sd,
    param_types  = fix$param_types
  )

  # Before correction, corrected_draws is NULL
  mat <- as.matrix(result)
  expect_true(is.matrix(mat))
  expect_identical(mat, result$original_draws)
})
