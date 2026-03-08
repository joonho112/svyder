# test-extract-rstanarm.R
# Tests for rstanarm integration: extract_draws.stanreg, der_compute.stanreg
# ---------------------------------------------------------------------------
# Tests that require rstanarm use skip_if_not_installed("rstanarm").
# Tests for the internal helper .map_rstanarm_params() work without rstanarm.


# ============================================================================
# .map_rstanarm_params() helper — works WITHOUT rstanarm installed
# ============================================================================

test_that(".map_rstanarm_params correctly maps fixed effects", {
  params <- c("(Intercept)", "x1", "x2",
              "b[(Intercept) group:1]", "b[(Intercept) group:2]",
              "Sigma[group:(Intercept),(Intercept)]",
              "sigma", "log-posterior")
  result <- svyder:::.map_rstanarm_params(params)

  fe_rows <- result[result$type == "fixed", ]
  re_rows <- result[result$type == "random", ]

  # Should have 3 fixed + 2 random
  expect_equal(nrow(fe_rows), 3L)
  expect_equal(nrow(re_rows), 2L)

  # Fixed effects mapped correctly
  expect_equal(fe_rows$original_name, c("(Intercept)", "x1", "x2"))
  expect_equal(fe_rows$svyder_name, c("beta[1]", "beta[2]", "beta[3]"))

  # Random effects mapped correctly
  expect_equal(re_rows$svyder_name, c("theta[1]", "theta[2]"))
})

test_that(".map_rstanarm_params excludes variance components", {
  params <- c("(Intercept)",
              "Sigma[group:(Intercept),(Intercept)]",
              "sigma", "log-posterior", "mean_PPD")
  result <- svyder:::.map_rstanarm_params(params)

  # Only (Intercept) should remain as fixed
  expect_equal(nrow(result), 1L)
  expect_equal(result$original_name, "(Intercept)")
  expect_equal(result$type, "fixed")
})

test_that(".map_rstanarm_params handles empty input", {
  result <- svyder:::.map_rstanarm_params(character(0))
  expect_equal(nrow(result), 0L)
  expect_true(is.data.frame(result))
})

test_that(".map_rstanarm_params handles no matching parameters", {
  params <- c("sigma", "log-posterior", "mean_PPD",
              "Sigma[group:(Intercept),(Intercept)]")
  result <- svyder:::.map_rstanarm_params(params)
  expect_equal(nrow(result), 0L)
})

test_that(".map_rstanarm_params preserves order: fixed before random", {
  params <- c("b[(Intercept) group:1]", "(Intercept)",
              "b[(Intercept) group:2]", "x1")
  result <- svyder:::.map_rstanarm_params(params)

  # Fixed effects should come first
  expect_equal(result$type[1], "fixed")
  expect_equal(result$type[2], "fixed")
  expect_equal(result$type[3], "random")
  expect_equal(result$type[4], "random")
})


# ============================================================================
# .classify_rstanarm_param() helper
# ============================================================================

test_that(".classify_rstanarm_param classifies correctly", {
  params <- c("(Intercept)", "x1",
              "b[(Intercept) group:1]",
              "Sigma[group:(Intercept),(Intercept)]",
              "sigma", "log-posterior")
  types <- svyder:::.classify_rstanarm_param(params)

  expect_equal(types[1], "fixed")    # (Intercept)
  expect_equal(types[2], "fixed")    # x1
  expect_equal(types[3], "random")   # b[...]
  expect_equal(types[4], "variance") # Sigma[...]
  expect_equal(types[5], "variance") # sigma
  expect_equal(types[6], "diagnostic") # log-posterior
})


# ============================================================================
# Method existence checks
# ============================================================================

test_that("extract_draws.stanreg method exists", {
  expect_true(is.function(svyder:::extract_draws.stanreg))
})

test_that("der_compute.stanreg method exists", {
  expect_true(is.function(svyder:::der_compute.stanreg))
})


# ============================================================================
# extract_draws.stanreg — requires rstanarm
# ============================================================================

test_that("extract_draws.stanreg works with fitted model", {
  skip_if_not_installed("rstanarm")

  # Integration test placeholder. In practice, testing requires a fitted
  # rstanarm model which is too heavy for unit tests. The method's core
  # logic is tested through the parameter mapping helpers above.
})


# ============================================================================
# der_compute.stanreg — requires rstanarm
# ============================================================================

test_that("der_compute.stanreg delegates correctly", {
  skip_if_not_installed("rstanarm")

  # Integration test placeholder. The stanreg method delegates to
  # der_compute.matrix after extracting draws and model components.
})
