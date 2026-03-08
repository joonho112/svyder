###############################################################################
# test-utils.R
# Tests for input validation and helper utilities (R/utils.R)
###############################################################################

# --- validate_inputs ---

test_that("validate_inputs passes with correct inputs", {
  N <- 20
  p <- 2
  J <- 4
  M <- 100
  y <- rbinom(N, 1, 0.5)
  X <- cbind(1, rnorm(N))
  group <- rep(1:J, each = N / J)
  weights <- runif(N, 0.5, 2)
  draws_beta <- matrix(rnorm(M * p), M, p)
  draws_theta <- matrix(rnorm(M * J), M, J)

  expect_true(.validate_inputs(y, X, group, weights, "binomial",
                               draws_beta, draws_theta))
})

test_that("validate_inputs catches dimension mismatch y vs X", {
  y <- 1:10
  X <- matrix(1, nrow = 11, ncol = 1)  # wrong N
  expect_error(
    .validate_inputs(y, X, rep(1, 10), rep(1, 10), "binomial",
                     matrix(1, 100, 1), matrix(1, 100, 1)),
    "mismatch"
  )
})

test_that("validate_inputs catches dimension mismatch y vs group", {
  N <- 10
  y <- rep(0, N)
  X <- matrix(1, N, 1)
  group <- rep(1, N + 1)  # wrong length
  expect_error(
    .validate_inputs(y, X, group, rep(1, N), "binomial",
                     matrix(1, 100, 1), matrix(1, 100, 1)),
    "mismatch"
  )
})

test_that("validate_inputs catches dimension mismatch y vs weights", {
  N <- 10
  y <- rep(0, N)
  X <- matrix(1, N, 1)
  expect_error(
    .validate_inputs(y, X, rep(1, N), rep(1, N + 1), "binomial",
                     matrix(1, 100, 1), matrix(1, 100, 1)),
    "mismatch"
  )
})

test_that("validate_inputs catches negative weights", {
  expect_error(
    .validate_inputs(rep(0, 10), matrix(1, 10, 1), rep(1, 10),
                     c(1, -1, rep(1, 8)), "binomial",
                     matrix(1, 100, 1), matrix(1, 100, 1)),
    "positive"
  )
})

test_that("validate_inputs catches zero weights", {
  expect_error(
    .validate_inputs(rep(0, 10), matrix(1, 10, 1), rep(1, 10),
                     c(0, rep(1, 9)), "binomial",
                     matrix(1, 100, 1), matrix(1, 100, 1)),
    "positive"
  )
})

test_that("validate_inputs catches invalid family", {
  expect_error(
    .validate_inputs(rep(0, 10), matrix(1, 10, 1), rep(1, 10),
                     rep(1, 10), "invalid_family",
                     matrix(1, 100, 1), matrix(1, 100, 1)),
    "family"
  )
})

test_that("validate_inputs catches non-sequential groups", {
  # Groups 1, 3 (missing 2)
  y <- rep(0, 10)
  X <- matrix(1, 10, 1)
  group <- c(rep(1, 5), rep(3, 5))
  expect_error(
    .validate_inputs(y, X, group, rep(1, 10), "binomial",
                     matrix(1, 100, 1), matrix(1, 100, 3)),
    "sequential"
  )
})

test_that("validate_inputs catches draws_beta col mismatch", {
  N <- 10
  p <- 2
  J <- 2
  M <- 50
  y <- rep(0, N)
  X <- matrix(1, N, p)
  group <- rep(1:J, each = N / J)
  # draws_beta has wrong number of columns
  expect_error(
    .validate_inputs(y, X, group, rep(1, N), "binomial",
                     matrix(1, M, p + 1), matrix(1, M, J)),
    "mismatch"
  )
})

test_that("validate_inputs catches draws_theta col mismatch", {
  N <- 10
  p <- 1
  J <- 2
  M <- 50
  y <- rep(0, N)
  X <- matrix(1, N, p)
  group <- rep(1:J, each = N / J)
  # draws_theta has wrong number of columns
  expect_error(
    .validate_inputs(y, X, group, rep(1, N), "binomial",
                     matrix(1, M, p), matrix(1, M, J + 1)),
    "mismatch"
  )
})

test_that("validate_inputs catches draws row mismatch", {
  N <- 10
  p <- 1
  J <- 2
  M <- 50
  y <- rep(0, N)
  X <- matrix(1, N, p)
  group <- rep(1:J, each = N / J)
  # Different number of rows in draws_beta vs draws_theta
  expect_error(
    .validate_inputs(y, X, group, rep(1, N), "binomial",
                     matrix(1, M, p), matrix(1, M + 10, J)),
    "mismatch"
  )
})

test_that("validate_inputs catches non-numeric y", {
  expect_error(
    .validate_inputs(letters[1:10], matrix(1, 10, 1), rep(1, 10),
                     rep(1, 10), "binomial",
                     matrix(1, 100, 1), matrix(1, 100, 1)),
    "numeric"
  )
})

test_that("validate_inputs catches non-matrix X", {
  expect_error(
    .validate_inputs(rep(0, 10), rep(1, 10), rep(1, 10),
                     rep(1, 10), "binomial",
                     matrix(1, 100, 1), matrix(1, 100, 1)),
    "matrix"
  )
})

test_that("validate_inputs accepts gaussian family", {
  N <- 20
  p <- 2
  J <- 4
  M <- 100
  y <- rnorm(N)
  X <- cbind(1, rnorm(N))
  group <- rep(1:J, each = N / J)
  weights <- runif(N, 0.5, 2)
  draws_beta <- matrix(rnorm(M * p), M, p)
  draws_theta <- matrix(rnorm(M * J), M, J)

  expect_true(.validate_inputs(y, X, group, weights, "gaussian",
                               draws_beta, draws_theta))
})


# --- validate_svyder ---

test_that("validate_svyder rejects non-svyder objects", {
  expect_error(.validate_svyder(list(a = 1)), "svyder")
})

test_that("validate_svyder detects missing components", {
  obj <- structure(list(der = 1), class = "svyder")
  expect_error(.validate_svyder(obj), "missing")
})

test_that("validate_svyder passes on complete object", {
  obj <- structure(
    list(
      der = 1,
      classification = data.frame(),
      V_sand = matrix(1),
      sigma_mcmc = matrix(1),
      diagnostics = list()
    ),
    class = "svyder"
  )
  expect_true(.validate_svyder(obj))
})


# --- compute_ci ---

test_that("compute_ci returns correct quantiles", {
  # T5 tolerance (statistical estimates)
  set.seed(1)
  draws <- matrix(rnorm(10000), ncol = 1)
  ci <- .compute_ci(draws, prob = 0.95)

  expect_equal(unname(ci[1, "lower"]), qnorm(0.025), tolerance = 0.05)
  expect_equal(unname(ci[1, "upper"]), qnorm(0.975), tolerance = 0.05)
})

test_that("compute_ci respects prob argument", {
  set.seed(42)
  draws <- matrix(rnorm(50000), ncol = 1)
  ci_90 <- .compute_ci(draws, prob = 0.90)
  ci_95 <- .compute_ci(draws, prob = 0.95)

  # 95% CI should be wider than 90% CI

  expect_true(ci_95[1, "upper"] - ci_95[1, "lower"] >
              ci_90[1, "upper"] - ci_90[1, "lower"])
})

test_that("compute_ci handles multiple parameters", {
  set.seed(123)
  draws <- matrix(rnorm(1000 * 3), ncol = 3)
  ci <- .compute_ci(draws, prob = 0.95)

  expect_equal(nrow(ci), 3)
  expect_equal(ncol(ci), 2)
  expect_equal(colnames(ci), c("lower", "upper"))
})

test_that("compute_ci handles vector input", {
  set.seed(1)
  draws <- rnorm(10000)
  ci <- .compute_ci(draws, prob = 0.95)

  expect_equal(nrow(ci), 1)
  expect_equal(ncol(ci), 2)
})


# --- format_pct ---

test_that("format_pct produces correct output", {
  expect_equal(.format_pct(0.019), "1.9%")
  expect_equal(.format_pct(0.5), "50.0%")
  expect_equal(.format_pct(1.0), "100.0%")
})

test_that("format_pct respects digits argument", {
  expect_equal(.format_pct(0.12345, digits = 2), "12.35%")
  expect_equal(.format_pct(0.12345, digits = 0), "12%")
})

test_that("format_pct handles zero and negative", {
  expect_equal(.format_pct(0), "0.0%")
  expect_equal(.format_pct(-0.05), "-5.0%")
})

test_that("format_pct vectorized", {
  result <- .format_pct(c(0.1, 0.5, 0.9))
  expect_equal(result, c("10.0%", "50.0%", "90.0%"))
  expect_length(result, 3)
})
