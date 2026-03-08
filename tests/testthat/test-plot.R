# test-plot.R
# Tests for plot.svyder(), autoplot.svyder(), and internal plot functions


# ============================================================================
# Helper: create a fully processed svyder object
# ============================================================================
.make_test_result <- function(correct = TRUE) {
  fix <- make_unbalanced_binomial()
  draws <- cbind(fix$draws_beta, fix$draws_theta)
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
  result
}


# ============================================================================
# plot.svyder profile
# ============================================================================

test_that("plot.svyder profile runs without error", {
  result <- .make_test_result()
  expect_silent(plot(result, type = "profile"))
})

test_that("plot.svyder profile returns ggplot when ggplot2 available", {
  skip_if_not_installed("ggplot2")
  result <- .make_test_result()
  p <- plot(result, type = "profile")
  expect_s3_class(p, "ggplot")
})

test_that("plot.svyder profile works without classification", {
  fix <- make_unbalanced_binomial()
  draws <- cbind(fix$draws_beta, fix$draws_theta)
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
  # No classification yet - should still work
  expect_silent(plot(result, type = "profile"))
})

test_that("plot.svyder profile respects custom tau", {
  skip_if_not_installed("ggplot2")
  result <- .make_test_result()
  p <- plot(result, type = "profile", tau = 2.0)
  expect_s3_class(p, "ggplot")
})


# ============================================================================
# plot.svyder decomposition
# ============================================================================

test_that("plot.svyder decomposition runs without error", {
  result <- .make_test_result()
  expect_silent(plot(result, type = "decomposition"))
})

test_that("plot.svyder decomposition returns ggplot", {
  skip_if_not_installed("ggplot2")
  result <- .make_test_result()
  p <- plot(result, type = "decomposition")
  expect_s3_class(p, "ggplot")
})


# ============================================================================
# plot.svyder comparison
# ============================================================================

test_that("plot.svyder comparison works after correction", {
  result <- .make_test_result(correct = TRUE)
  expect_silent(plot(result, type = "comparison"))
})

test_that("plot.svyder comparison returns ggplot", {
  skip_if_not_installed("ggplot2")
  result <- .make_test_result(correct = TRUE)
  p <- plot(result, type = "comparison")
  expect_s3_class(p, "ggplot")
})

test_that("plot.svyder comparison errors without correction", {
  result <- .make_test_result(correct = FALSE)
  expect_error(plot(result, type = "comparison"), "corrected draws")
})

test_that("plot.svyder comparison works with explicit params", {
  skip_if_not_installed("ggplot2")
  result <- .make_test_result(correct = TRUE)
  p <- plot(result, type = "comparison", params = result$params[1:2])
  expect_s3_class(p, "ggplot")
})

test_that("plot.svyder comparison errors with invalid params", {
  result <- .make_test_result(correct = TRUE)
  expect_error(plot(result, type = "comparison", params = "nonexistent"),
               "No matching parameters")
})


# ============================================================================
# Base R fallback functions exist
# ============================================================================

test_that("base R fallback functions exist", {
  expect_true(is.function(svyder:::.plot_profile_base))
  expect_true(is.function(svyder:::.plot_decomposition_base))
  expect_true(is.function(svyder:::.plot_comparison_base))
})

test_that("base R profile plot returns invisible svyder", {
  result <- .make_test_result()
  ret <- svyder:::.plot_profile_base(result)
  expect_s3_class(ret, "svyder")
})

test_that("base R decomposition plot returns invisible svyder", {
  result <- .make_test_result()
  ret <- svyder:::.plot_decomposition_base(result)
  expect_s3_class(ret, "svyder")
})

test_that("base R comparison plot returns invisible svyder", {
  result <- .make_test_result(correct = TRUE)
  ret <- svyder:::.plot_comparison_base(result)
  expect_s3_class(ret, "svyder")
})


# ============================================================================
# autoplot.svyder
# ============================================================================

test_that("autoplot.svyder profile works", {
  skip_if_not_installed("ggplot2")
  result <- .make_test_result()
  p <- ggplot2::autoplot(result, type = "profile")
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.svyder decomposition works", {
  skip_if_not_installed("ggplot2")
  result <- .make_test_result()
  p <- ggplot2::autoplot(result, type = "decomposition")
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.svyder comparison works", {
  skip_if_not_installed("ggplot2")
  result <- .make_test_result(correct = TRUE)
  p <- ggplot2::autoplot(result, type = "comparison")
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.svyder defaults to profile", {
  skip_if_not_installed("ggplot2")
  result <- .make_test_result()
  p <- ggplot2::autoplot(result)
  expect_s3_class(p, "ggplot")
})
