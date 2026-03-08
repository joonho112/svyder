###############################################################################
# test-shrinkage.R
# Tests for per-group design effect and shrinkage computations (R/shrinkage.R)
###############################################################################

# --- compute_deff_j ---

test_that("equal weights give DEFF = 1.0", {
  # T1 tolerance (exact identity)
  group <- rep(1:5, each = 20)
  w <- rep(1, 100)
  deff <- .compute_deff_j(group, w)
  expect_equal(unname(deff), rep(1.0, 5), tolerance = 1e-12)
})

test_that("DEFF names are correct", {
  group <- rep(1:3, each = 10)
  w <- rep(1, 30)
  deff <- .compute_deff_j(group, w)
  expect_equal(names(deff), c("group_1", "group_2", "group_3"))
})

test_that("unequal weights give DEFF > 1", {
  group <- rep(1, 10)
  w <- c(rep(1, 5), rep(5, 5))
  deff <- .compute_deff_j(group, w)
  expect_true(deff > 1.0)
})

test_that("DEFF with known analytic result", {
  # For 2 observations with weights w1 and w2:
  # DEFF = 2 * (w1^2 + w2^2) / (w1 + w2)^2
  group <- rep(1, 2)
  w <- c(1, 3)
  deff <- .compute_deff_j(group, w)
  expected <- 2 * (1 + 9) / (4)^2  # 2 * 10 / 16 = 1.25

  expect_equal(unname(deff), expected, tolerance = 1e-12)
})

test_that("DEFF computed independently per group", {
  group <- c(rep(1, 10), rep(2, 10))
  w <- c(rep(1, 10), c(rep(1, 5), rep(5, 5)))
  deff <- .compute_deff_j(group, w)
  # Group 1: equal weights -> DEFF = 1
  expect_equal(unname(deff[1]), 1.0, tolerance = 1e-12)
  # Group 2: unequal weights -> DEFF > 1
  expect_true(deff[2] > 1.0)
})


# --- compute_B_j ---

test_that("B_j in [0, 1]", {
  group <- rep(1:5, each = 20)
  w <- rep(1, 100)
  wt <- rep(0.25, 100)  # mu=0.5
  sigma2 <- 0.25
  B <- .compute_B_j(group, w, wt, sigma2)
  expect_true(all(B >= 0 & B <= 1))
})

test_that("B_j names are correct", {
  group <- rep(1:3, each = 10)
  w <- rep(1, 30)
  wt <- rep(0.25, 30)
  sigma2 <- 1.0
  B <- .compute_B_j(group, w, wt, sigma2)
  expect_equal(names(B), c("group_1", "group_2", "group_3"))
})

test_that("B_j increases with sigma2", {
  group <- rep(1:5, each = 20)
  w <- rep(1, 100)
  wt <- rep(0.25, 100)
  B_small <- .compute_B_j(group, w, wt, sigma2 = 0.01)
  B_large <- .compute_B_j(group, w, wt, sigma2 = 10.0)
  expect_true(all(B_large > B_small))
})

test_that("B_j approaches 1 for very large sigma2", {
  group <- rep(1:3, each = 10)
  w <- rep(1, 30)
  wt <- rep(0.25, 30)
  B <- .compute_B_j(group, w, wt, sigma2 = 1e6)
  expect_true(all(B > 0.999))
})

test_that("B_j approaches 0 for very small sigma2", {
  group <- rep(1:3, each = 10)
  w <- rep(1, 30)
  wt <- rep(0.25, 30)
  B <- .compute_B_j(group, w, wt, sigma2 = 1e-10)
  expect_true(all(B < 0.001))
})

test_that("B_j analytic formula: sigma2 / (sigma2 + V_tilde_j)", {
  # T3 tolerance (analytical formula)
  group <- rep(1, 10)
  w <- rep(2, 10)
  wt <- rep(0.25, 10)
  sigma2 <- 0.5

  B <- .compute_B_j(group, w, wt, sigma2)

  # info_j = sum(w * wt) = 10 * 2 * 0.25 = 5
  # V_tilde_j = 1/5 = 0.2
  # B_j = 0.5 / (0.5 + 0.2) = 0.5/0.7
  expected <- 0.5 / 0.7
  expect_equal(unname(B), expected, tolerance = 1e-6)
})


# --- compute_kappa_j ---

test_that("kappa(J) in [0, 1]", {
  B <- c(0.2, 0.5, 0.8, 0.95)
  J <- 10
  kappa <- .compute_kappa_j(B, J)
  expect_true(all(kappa >= 0 & kappa <= 1))
})

test_that("kappa = 0 when B = 1", {
  expect_equal(.compute_kappa_j(1.0, 10), 0.0)
})

test_that("kappa approaches (J-1)/J as B -> 0", {
  J <- 100
  kappa <- .compute_kappa_j(0.001, J)
  expect_equal(kappa, (J - 1) / J, tolerance = 0.01)
})

test_that("kappa at B = 0 is exactly (J-1)/J", {
  # T1 tolerance (exact identity)
  J <- 50
  kappa <- .compute_kappa_j(0.0, J)
  expect_equal(kappa, (J - 1) / J, tolerance = 1e-12)
})

test_that("kappa is monotonically decreasing in B", {
  B <- seq(0, 1, by = 0.1)
  J <- 20
  kappa <- .compute_kappa_j(B, J)
  # Each subsequent value should be <= previous
  diffs <- diff(kappa)
  expect_true(all(diffs <= 0))
})

test_that("kappa vectorized over B_j", {
  B <- c(0.0, 0.5, 1.0)
  J <- 10
  kappa <- .compute_kappa_j(B, J)
  expect_length(kappa, 3)
  expect_equal(kappa[1], (J - 1) / J, tolerance = 1e-12)
  expect_equal(kappa[3], 0.0, tolerance = 1e-12)
})
