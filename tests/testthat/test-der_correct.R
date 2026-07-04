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
    cluster      = fix$psu,
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
    cluster      = fix$psu,
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
    cluster      = fix$psu,
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
    cluster      = fix$psu,
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
    cluster      = fix$psu,
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
    cluster      = fix$psu,
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
    cluster      = fix$psu,
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


# ============================================================================
# Correction methods: block_cholesky vs marginal
# ============================================================================

# Helper: classified object with at least two flagged parameters
.classified_two_flagged <- function() {
  fix <- make_unbalanced_binomial()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)
  result <- der_compute(
    draws_all,
    y            = fix$y,
    X            = fix$X,
    group        = fix$group,
    weights      = fix$w,
    cluster      = fix$psu,
    family       = "binomial",
    sigma_theta  = fix$sigma_theta_hat,
    beta_prior_sd = fix$beta_prior_sd,
    param_types  = fix$param_types
  )
  # Force the two fixed effects to be flagged
  result$der[1] <- 2.0
  result$der[2] <- 1.8
  der_classify(result, tau = 1.2, verbose = FALSE)
}

test_that("block_cholesky matches the full sandwich block covariance", {
  result <- der_correct(.classified_two_flagged(), method = "block_cholesky")

  F_idx <- which(result$classification$flagged)
  expect_gte(length(F_idx), 2L)

  cov_corrected <- cov(result$corrected_draws[, F_idx, drop = FALSE])
  V_F           <- result$V_sand[F_idx, F_idx, drop = FALSE]
  expect_equal(unname(cov_corrected), unname(V_F), tolerance = 1e-8)

  # Means preserved
  expect_equal(colMeans(result$corrected_draws),
               colMeans(result$original_draws), tolerance = 1e-10)

  # Unflagged columns bitwise identical
  for (i in setdiff(seq_along(result$der), F_idx)) {
    expect_identical(result$corrected_draws[, i], result$original_draws[, i])
  }

  expect_equal(result$correction_method, "block_cholesky")
})

test_that("marginal method matches per-parameter sqrt(DER) scaling", {
  classified <- .classified_two_flagged()
  result <- der_correct(classified, method = "marginal")

  F_idx     <- which(result$classification$flagged)
  point_est <- colMeans(result$original_draws)

  for (i in F_idx) {
    sf <- sqrt(diag(result$V_sand)[i] / diag(result$sigma_mcmc)[i])
    expected_col <- point_est[i] +
      sf * (result$original_draws[, i] - point_est[i])
    expect_equal(result$corrected_draws[, i], expected_col, tolerance = 1e-12)
    # Marginal variance matches the target for each flagged parameter
    expect_equal(cov(result$corrected_draws[, i, drop = FALSE])[1, 1],
                 diag(result$V_sand)[i], tolerance = 1e-8)
  }

  expect_equal(result$correction_method, "marginal")
})

test_that("block and marginal differ off-diagonally with >= 2 flagged", {
  classified <- .classified_two_flagged()
  res_block    <- der_correct(classified, method = "block_cholesky")
  res_marginal <- der_correct(classified, method = "marginal")

  F_idx <- which(classified$classification$flagged)

  # Marginal SDs agree between methods...
  expect_equal(apply(res_block$corrected_draws[, F_idx], 2, sd),
               apply(res_marginal$corrected_draws[, F_idx], 2, sd),
               tolerance = 1e-8)
  # ...but the joint draws differ (marginal does not match the covariance)
  expect_false(isTRUE(all.equal(res_block$corrected_draws[, F_idx],
                                res_marginal$corrected_draws[, F_idx],
                                tolerance = 1e-10)))

  # Marginal leaves the flagged-block correlation at its MCMC value,
  # so its cross-covariance need not equal V_sand off-diagonals
  cov_marg <- cov(res_marginal$corrected_draws[, F_idx])
  V_F      <- classified$V_sand[F_idx, F_idx]
  expect_false(isTRUE(all.equal(unname(cov_marg), unname(V_F),
                                tolerance = 1e-6)))
})

test_that("methods coincide when exactly one parameter is flagged", {
  fix <- make_unbalanced_binomial()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)
  result <- der_compute(
    draws_all,
    y = fix$y, X = fix$X, group = fix$group,
    weights = fix$w, cluster = fix$psu,
    family = "binomial", sigma_theta = fix$sigma_theta_hat,
    beta_prior_sd = fix$beta_prior_sd, param_types = fix$param_types
  )
  result$der[] <- 0.5
  result$der[1] <- 3.0   # single flagged parameter
  result <- der_classify(result, tau = 1.2, verbose = FALSE)

  res_block    <- der_correct(result, method = "block_cholesky")
  res_marginal <- der_correct(result, method = "marginal")

  expect_equal(res_block$corrected_draws, res_marginal$corrected_draws,
               tolerance = 1e-12)
})

test_that("method = 'cholesky' is retired with an informative error", {
  classified <- .classified_two_flagged()
  expect_error(der_correct(classified, method = "cholesky"),
               "retired")
})
