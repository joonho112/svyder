# test-extract.R
# Tests for the core extract_draws() and extract_design() generics
# ---------------------------------------------------------------------------
# These tests do NOT require Stan, brms, or other optional backends.


# ============================================================================
# extract_draws.matrix
# ============================================================================

test_that("extract_draws.matrix is identity", {
  m <- matrix(1:20, nrow = 5)
  result <- extract_draws(m)
  expect_type(result, "list")
  expect_identical(result$draws, m)
})

test_that("extract_draws.matrix preserves column names", {
  m <- matrix(rnorm(30), nrow = 10, ncol = 3)
  colnames(m) <- c("beta[1]", "theta[1]", "theta[2]")
  result <- extract_draws(m)
  expect_identical(colnames(result$draws), c("beta[1]", "theta[1]", "theta[2]"))
})

test_that("extract_draws.matrix works with single column", {
  m <- matrix(rnorm(10), ncol = 1)
  result <- extract_draws(m)
  expect_equal(ncol(result$draws), 1L)
  expect_equal(nrow(result$draws), 10L)
})

test_that("extract_draws.matrix works with single row", {
  m <- matrix(1:5, nrow = 1)
  result <- extract_draws(m)
  expect_equal(nrow(result$draws), 1L)
  expect_equal(ncol(result$draws), 5L)
})


# ============================================================================
# extract_draws.default
# ============================================================================

test_that("extract_draws.default errors informatively on list", {
  expect_error(
    extract_draws(list(a = 1, b = 2)),
    "does not know how to handle class 'list'"
  )
})

test_that("extract_draws.default errors informatively on data.frame", {
  expect_error(
    extract_draws(data.frame(x = 1:3)),
    "does not know how to handle class 'data.frame'"
  )
})

test_that("extract_draws.default errors on character", {
  expect_error(
    extract_draws("not_a_model"),
    "does not know how to handle class 'character'"
  )
})

test_that("extract_draws.default mentions supported classes", {
  expect_error(
    extract_draws(42L),
    "Supported: matrix, brmsfit, CmdStanMCMC, stanreg, draws_matrix, draws_df"
  )
})


# ============================================================================
# extract_design.default
# ============================================================================

test_that("extract_design.default errors informatively on list", {
  expect_error(
    extract_design(list()),
    "does not support class 'list'"
  )
})

test_that("extract_design.default errors informatively on data.frame", {
  expect_error(
    extract_design(data.frame(x = 1)),
    "does not support class 'data.frame'"
  )
})

test_that("extract_design.default mentions supported classes", {
  expect_error(
    extract_design(42),
    "Supported: survey.design2"
  )
})


# ============================================================================
# extract_design.survey.design2
# ============================================================================

test_that("extract_design works with survey.design2 object", {
  skip_if_not_installed("survey")

  # Create a simple survey design
  dat <- data.frame(
    y       = rnorm(20),
    weights = rep(c(1.5, 0.8), each = 10),
    psu     = rep(1:4, each = 5),
    strata  = rep(c("A", "B"), each = 10)
  )

  design <- survey::svydesign(
    ids     = ~psu,
    strata  = ~strata,
    weights = ~weights,
    data    = dat
  )

  result <- extract_design(design)

  expect_type(result, "list")
  expect_true("weights" %in% names(result))
  expect_true("cluster" %in% names(result))
  expect_true("strata" %in% names(result))
  expect_equal(length(result$weights), 20L)
})

test_that("extract_design preserves survey weights correctly", {
  skip_if_not_installed("survey")

  w <- c(1.2, 0.8, 1.5, 0.6, 1.1)
  dat <- data.frame(
    y       = rnorm(5),
    wt      = w,
    cluster = 1:5
  )

  design <- survey::svydesign(
    ids     = ~cluster,
    weights = ~wt,
    data    = dat
  )

  result <- extract_design(design)
  expect_equal(as.numeric(result$weights), w)
})
