# test-tidy.R
# Tests for tidy.svyder() and glance.svyder()


# ============================================================================
# Helper
# ============================================================================
.make_tidy_result <- function(classify = TRUE, correct = TRUE) {
  fix <- make_unbalanced_binomial()
  draws <- cbind(fix$draws_beta, fix$draws_theta)

  if (classify) {
    result <- der_diagnose(
      draws,
      y            = fix$y,
      X            = fix$X,
      group        = fix$group,
      weights      = fix$w,
      family       = "binomial",
      sigma_theta  = fix$sigma_theta_hat,
      param_types  = fix$param_types,
      tau          = 1.2,
      correct      = correct
    )
  } else {
    result <- der_compute(
      draws,
      y            = fix$y,
      X            = fix$X,
      group        = fix$group,
      weights      = fix$w,
      family       = "binomial",
      sigma_theta  = fix$sigma_theta_hat,
      param_types  = fix$param_types
    )
  }
  result
}


# ============================================================================
# tidy.svyder
# ============================================================================

test_that("tidy returns a data.frame", {
  result <- .make_tidy_result()
  td <- tidy.svyder(result)
  expect_s3_class(td, "data.frame")
})

test_that("tidy returns one row per parameter", {
  result <- .make_tidy_result()
  td <- tidy.svyder(result)
  expect_equal(nrow(td), length(result$der))
})

test_that("tidy has correct columns with classification", {
  result <- .make_tidy_result()
  td <- tidy.svyder(result)
  expected_cols <- c("term", "estimate", "std.error", "der",
                     "tier", "action", "flagged", "scale_factor")
  expect_true(all(expected_cols %in% names(td)))
})

test_that("tidy has correct columns without classification", {
  result <- .make_tidy_result(classify = FALSE)
  td <- tidy.svyder(result)
  # Must always include term and der
  expect_true("term" %in% names(td))
  expect_true("der" %in% names(td))
  expect_true("estimate" %in% names(td))
  expect_true("std.error" %in% names(td))
})

test_that("tidy term matches param names", {
  result <- .make_tidy_result()
  td <- tidy.svyder(result)
  expect_equal(td$term, result$params)
})

test_that("tidy der values match svyder$der", {
  result <- .make_tidy_result()
  td <- tidy.svyder(result)
  expect_equal(td$der, as.numeric(result$der))
})

test_that("tidy estimate is posterior mean from original draws", {
  result <- .make_tidy_result()
  td <- tidy.svyder(result)
  expected_mean <- unname(colMeans(result$original_draws))
  expect_equal(td$estimate, expected_mean, tolerance = 1e-10)
})

test_that("tidy std.error is posterior SD", {
  result <- .make_tidy_result()
  td <- tidy.svyder(result)
  expected_sd <- unname(apply(result$original_draws, 2, sd))
  expect_equal(td$std.error, expected_sd, tolerance = 1e-10)
})

test_that("tidy flagged column has correct type", {
  result <- .make_tidy_result()
  td <- tidy.svyder(result)
  expect_type(td$flagged, "logical")
})

test_that("tidy scale_factor reflects corrections", {
  result <- .make_tidy_result(correct = TRUE)
  td <- tidy.svyder(result)
  # Unflagged parameters should have scale_factor = 1.0
  unflagged <- which(!td$flagged)
  if (length(unflagged) > 0) {
    expect_true(all(td$scale_factor[unflagged] == 1.0))
  }
})

test_that("tidy without classification has NA tier/action/flagged", {
  result <- .make_tidy_result(classify = FALSE)
  td <- tidy.svyder(result)
  expect_true(all(is.na(td$tier)))
  expect_true(all(is.na(td$action)))
  expect_true(all(is.na(td$flagged)))
})

test_that("generics::tidy dispatches correctly when generics loaded", {
  skip_if_not_installed("generics")
  result <- .make_tidy_result()
  td <- generics::tidy(result)
  expect_s3_class(td, "data.frame")
  expect_true("term" %in% names(td))
  expect_true("der" %in% names(td))
})


# ============================================================================
# glance.svyder
# ============================================================================

test_that("glance returns a one-row data.frame", {
  result <- .make_tidy_result()
  gl <- glance.svyder(result)
  expect_s3_class(gl, "data.frame")
  expect_equal(nrow(gl), 1)
})

test_that("glance has correct columns", {
  result <- .make_tidy_result()
  gl <- glance.svyder(result)
  expected_cols <- c("n_params", "n_flagged", "tau", "family",
                     "n_obs", "n_groups", "mean_deff", "mean_B",
                     "der_min", "der_max")
  expect_true(all(expected_cols %in% names(gl)))
})

test_that("glance n_params matches parameter count", {
  result <- .make_tidy_result()
  gl <- glance.svyder(result)
  expect_equal(gl$n_params, length(result$der))
})

test_that("glance n_flagged matches classification", {
  result <- .make_tidy_result()
  gl <- glance.svyder(result)
  expected <- sum(result$classification$flagged)
  expect_equal(gl$n_flagged, expected)
})

test_that("glance tau matches stored tau", {
  result <- .make_tidy_result()
  gl <- glance.svyder(result)
  expect_equal(gl$tau, result$tau)
})

test_that("glance family matches stored family", {
  result <- .make_tidy_result()
  gl <- glance.svyder(result)
  expect_equal(gl$family, result$family)
})

test_that("glance der_min and der_max are correct", {
  result <- .make_tidy_result()
  gl <- glance.svyder(result)
  expect_equal(gl$der_min, min(result$der))
  expect_equal(gl$der_max, max(result$der))
})

test_that("glance mean_deff is mean of deff_j", {
  result <- .make_tidy_result()
  gl <- glance.svyder(result)
  expect_equal(gl$mean_deff, mean(result$deff_j))
})

test_that("glance mean_B is mean of B_j", {
  result <- .make_tidy_result()
  gl <- glance.svyder(result)
  expect_equal(gl$mean_B, mean(result$B_j))
})

test_that("glance without classification has NA n_flagged", {
  result <- .make_tidy_result(classify = FALSE)
  gl <- glance.svyder(result)
  expect_true(is.na(gl$n_flagged))
})

test_that("generics::glance dispatches correctly when generics loaded", {
  skip_if_not_installed("generics")
  result <- .make_tidy_result()
  gl <- generics::glance(result)
  expect_s3_class(gl, "data.frame")
  expect_equal(nrow(gl), 1)
  expect_true("n_params" %in% names(gl))
})
