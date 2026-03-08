# test-extract-brms.R
# Tests for brms integration: extract_draws.brmsfit, der_compute.brmsfit
# ---------------------------------------------------------------------------
# Tests that require brms use skip_if_not_installed("brms").
# Tests for the internal helper .map_brms_params() work without brms.


# ============================================================================
# .map_brms_params() helper — works WITHOUT brms installed
# ============================================================================

test_that(".map_brms_params correctly maps fixed effects", {
  params <- c("b_Intercept", "b_x1", "b_x2",
              "sd_group__Intercept", "sigma", "lp__")
  result <- svyder:::.map_brms_params(params)

  # Should only include b_ parameters (not sd_, sigma, lp__)
  expect_equal(nrow(result), 3L)
  expect_equal(result$original_name, c("b_Intercept", "b_x1", "b_x2"))
  expect_equal(result$svyder_name, c("beta[1]", "beta[2]", "beta[3]"))
  expect_true(all(result$type == "fixed"))
})

test_that(".map_brms_params correctly maps random effects", {
  params <- c("b_Intercept",
              "r_group[1,Intercept]", "r_group[2,Intercept]",
              "r_group[3,Intercept]",
              "sd_group__Intercept", "lp__")
  result <- svyder:::.map_brms_params(params)

  # Should have 1 fixed + 3 random

  expect_equal(nrow(result), 4L)

  fe_rows <- result[result$type == "fixed", ]
  re_rows <- result[result$type == "random", ]

  expect_equal(nrow(fe_rows), 1L)
  expect_equal(nrow(re_rows), 3L)

  expect_equal(fe_rows$svyder_name, "beta[1]")
  expect_equal(re_rows$svyder_name, c("theta[1]", "theta[2]", "theta[3]"))
})

test_that(".map_brms_params excludes variance components", {
  params <- c("b_Intercept", "b_x1",
              "sd_group__Intercept", "sd_group__x1",
              "cor_group__Intercept__x1",
              "sigma", "lp__", "lprior")
  result <- svyder:::.map_brms_params(params)

  # Only b_ parameters should be included
  expect_equal(nrow(result), 2L)
  expect_true(all(result$type == "fixed"))
  expect_equal(result$original_name, c("b_Intercept", "b_x1"))
})

test_that(".map_brms_params excludes distributional parameters", {
  params <- c("b_Intercept", "b_x1",
              "b_sigma_Intercept", "b_sigma_x1",
              "sd_group__Intercept")
  result <- svyder:::.map_brms_params(params)

  # b_sigma_ should be excluded
  expect_equal(nrow(result), 2L)
  expect_equal(result$original_name, c("b_Intercept", "b_x1"))
})

test_that(".map_brms_params handles empty input", {
  result <- svyder:::.map_brms_params(character(0))
  expect_equal(nrow(result), 0L)
  expect_true(is.data.frame(result))
})

test_that(".map_brms_params handles no matching parameters", {
  params <- c("sigma", "lp__", "sd_group__Intercept")
  result <- svyder:::.map_brms_params(params)
  expect_equal(nrow(result), 0L)
})

test_that(".map_brms_params preserves order: fixed before random", {
  params <- c("r_group[1,Intercept]", "b_Intercept",
              "r_group[2,Intercept]", "b_x1",
              "r_group[3,Intercept]")
  result <- svyder:::.map_brms_params(params)

  # Fixed effects should come first
  expect_equal(result$type[1], "fixed")
  expect_equal(result$type[2], "fixed")
  expect_equal(result$type[3], "random")
  expect_equal(result$type[4], "random")
  expect_equal(result$type[5], "random")
})

test_that(".map_brms_params handles complex random effect names", {
  params <- c("b_Intercept",
              "r_school[SchoolA,Intercept]",
              "r_school[SchoolB,Intercept]",
              "r_school[School With Spaces,Intercept]")
  result <- svyder:::.map_brms_params(params)

  re_rows <- result[result$type == "random", ]
  expect_equal(nrow(re_rows), 3L)
  expect_equal(re_rows$svyder_name, c("theta[1]", "theta[2]", "theta[3]"))
})


# ============================================================================
# extract_draws.brmsfit — requires brms
# ============================================================================

test_that("extract_draws.brmsfit requires brms", {
  skip_if_not_installed("brms")
  # If brms is installed, this test verifies the method exists

  expect_true(is.function(svyder:::extract_draws.brmsfit))
})


# ============================================================================
# der_compute.brmsfit — requires brms (integration test)
# ============================================================================

test_that("der_compute.brmsfit method exists", {
  # This test verifies the method is registered regardless of brms
  expect_true(is.function(svyder:::der_compute.brmsfit))
})
