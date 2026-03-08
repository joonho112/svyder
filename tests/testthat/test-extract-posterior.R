# test-extract-posterior.R
# Tests for posterior package integration: draws_matrix, draws_df, etc.
# ---------------------------------------------------------------------------
# All tests require the posterior package.


# ============================================================================
# extract_draws.draws_matrix
# ============================================================================

test_that("extract_draws.draws_matrix converts to plain matrix", {
  skip_if_not_installed("posterior")

  # Create a draws_matrix object
  m <- matrix(rnorm(30), nrow = 10, ncol = 3)
  colnames(m) <- c("alpha", "beta", "sigma")
  dm <- posterior::as_draws_matrix(m)

  result <- extract_draws(dm)

  expect_type(result, "list")
  expect_true(is.matrix(result$draws))
  expect_equal(ncol(result$draws), 3L)
  expect_equal(nrow(result$draws), 10L)
  expect_true(all(c("alpha", "beta", "sigma") %in% colnames(result$draws)))
})

test_that("extract_draws.draws_matrix removes metadata columns", {
  skip_if_not_installed("posterior")

  m <- matrix(rnorm(20), nrow = 10, ncol = 2)
  colnames(m) <- c("mu", "tau")
  dm <- posterior::as_draws_matrix(m)

  result <- extract_draws(dm)

  # .chain, .iteration, .draw should NOT be in the result
  expect_false(".chain" %in% colnames(result$draws))
  expect_false(".iteration" %in% colnames(result$draws))
  expect_false(".draw" %in% colnames(result$draws))
})


# ============================================================================
# extract_draws.draws_df
# ============================================================================

test_that("extract_draws.draws_df converts to plain matrix", {
  skip_if_not_installed("posterior")

  m <- matrix(rnorm(30), nrow = 10, ncol = 3)
  colnames(m) <- c("alpha", "beta", "sigma")
  ddf <- posterior::as_draws_df(m)

  result <- extract_draws(ddf)

  expect_type(result, "list")
  expect_true(is.matrix(result$draws))
  expect_equal(ncol(result$draws), 3L)
  expect_true(all(c("alpha", "beta", "sigma") %in% colnames(result$draws)))
})

test_that("extract_draws.draws_df removes metadata columns", {
  skip_if_not_installed("posterior")

  m <- matrix(rnorm(20), nrow = 10, ncol = 2)
  colnames(m) <- c("mu", "tau")
  ddf <- posterior::as_draws_df(m)

  result <- extract_draws(ddf)

  expect_false(".chain" %in% colnames(result$draws))
  expect_false(".iteration" %in% colnames(result$draws))
  expect_false(".draw" %in% colnames(result$draws))
})


# ============================================================================
# extract_draws.draws_array
# ============================================================================

test_that("extract_draws.draws_array merges chains", {
  skip_if_not_installed("posterior")

  # Create a draws_array with 2 chains x 5 iterations x 3 parameters
  arr <- array(rnorm(30), dim = c(5, 2, 3),
               dimnames = list(NULL, NULL, c("mu", "sigma", "tau")))
  da <- posterior::as_draws_array(arr)

  result <- extract_draws(da)

  expect_true(is.matrix(result$draws))
  # 2 chains x 5 iterations = 10 total draws
  expect_equal(nrow(result$draws), 10L)
  expect_equal(ncol(result$draws), 3L)
})


# ============================================================================
# extract_draws.draws_list
# ============================================================================

test_that("extract_draws.draws_list converts to matrix", {
  skip_if_not_installed("posterior")

  # Create a draws_list
  dl <- posterior::as_draws_list(list(
    list(mu = rnorm(5), sigma = abs(rnorm(5))),
    list(mu = rnorm(5), sigma = abs(rnorm(5)))
  ))

  result <- extract_draws(dl)

  expect_true(is.matrix(result$draws))
  # 2 chains x 5 iterations = 10 total draws
  expect_equal(nrow(result$draws), 10L)
  expect_equal(ncol(result$draws), 2L)
  expect_true(all(c("mu", "sigma") %in% colnames(result$draws)))
})


# ============================================================================
# extract_draws.draws_rvars
# ============================================================================

test_that("extract_draws.draws_rvars converts to matrix", {
  skip_if_not_installed("posterior")

  # Create draws_rvars from a draws_matrix
  m <- matrix(rnorm(20), nrow = 10, ncol = 2)
  colnames(m) <- c("mu", "sigma")
  drv <- posterior::as_draws_rvars(m)

  result <- extract_draws(drv)

  expect_true(is.matrix(result$draws))
  expect_equal(ncol(result$draws), 2L)
  expect_true(all(c("mu", "sigma") %in% colnames(result$draws)))
})


# ============================================================================
# Method existence checks
# ============================================================================

test_that("all posterior extract_draws methods exist", {
  expect_true(is.function(svyder:::extract_draws.draws_matrix))
  expect_true(is.function(svyder:::extract_draws.draws_df))
  expect_true(is.function(svyder:::extract_draws.draws_array))
  expect_true(is.function(svyder:::extract_draws.draws_list))
  expect_true(is.function(svyder:::extract_draws.draws_rvars))
})
