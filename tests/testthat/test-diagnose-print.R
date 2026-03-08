# test-diagnose-print.R
# Tests for der_diagnose(), print.svyder(), summary.svyder()


# ============================================================================
# der_diagnose tests
# ============================================================================

test_that("der_diagnose matches stepwise pipeline", {
  fix <- make_unbalanced_binomial()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)

  # Stepwise
  step_result <- der_compute(
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
  step_result <- der_classify(step_result, tau = 1.2, verbose = FALSE)
  step_result <- der_correct(step_result)

  # All-in-one
  diag_result <- der_diagnose(
    draws_all,
    y            = fix$y,
    X            = fix$X,
    group        = fix$group,
    weights      = fix$w,
    psu          = fix$psu,
    family       = "binomial",
    sigma_theta  = fix$sigma_theta_hat,
    beta_prior_sd = fix$beta_prior_sd,
    param_types  = fix$param_types,
    tau          = 1.2,
    correct      = TRUE
  )

  # DER values should be identical
  expect_equal(as.numeric(diag_result$der), as.numeric(step_result$der))

  # Classification should be identical
  expect_equal(diag_result$classification$flagged,
               step_result$classification$flagged)
  expect_equal(diag_result$classification$tier,
               step_result$classification$tier)
  expect_equal(diag_result$classification$action,
               step_result$classification$action)

  # Scale factors should be identical
  expect_equal(diag_result$scale_factors, step_result$scale_factors)
})

test_that("der_diagnose without correction", {
  fix <- make_balanced_gaussian()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)

  result <- der_diagnose(
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
    param_types  = fix$param_types,
    tau          = 1.2,
    correct      = FALSE
  )

  expect_s3_class(result, "svyder")
  # Should have classification but no corrected draws
  expect_true("flagged" %in% names(result$classification))
  expect_null(result$corrected_draws)
})


# ============================================================================
# print.svyder tests
# ============================================================================

test_that("print returns invisible svyder", {
  fix <- make_balanced_gaussian()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)

  result <- der_diagnose(
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
    param_types  = fix$param_types,
    tau          = 1.2
  )

  out <- withr::with_output_sink(tempfile(), print(result))
  expect_s3_class(out, "svyder")
})

test_that("print output contains 'svyder diagnostic'", {
  fix <- make_balanced_gaussian()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)

  result <- der_diagnose(
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
    param_types  = fix$param_types,
    tau          = 1.2
  )

  output <- capture.output(print(result))
  expect_true(any(grepl("svyder diagnostic", output)))
})

test_that("print works on unclassified svyder object", {
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

  output <- capture.output(print(result))
  expect_true(any(grepl("svyder diagnostic", output)))
  expect_true(any(grepl("not yet classified", output)))
})


# ============================================================================
# summary.svyder tests
# ============================================================================

test_that("summary returns data.frame", {
  fix <- make_unbalanced_binomial()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)

  result <- der_diagnose(
    draws_all,
    y            = fix$y,
    X            = fix$X,
    group        = fix$group,
    weights      = fix$w,
    psu          = fix$psu,
    family       = "binomial",
    sigma_theta  = fix$sigma_theta_hat,
    beta_prior_sd = fix$beta_prior_sd,
    param_types  = fix$param_types,
    tau          = 1.2
  )

  out <- withr::with_output_sink(tempfile(), summary(result))
  expect_true(is.data.frame(out))
})

test_that("summary has correct columns", {
  fix <- make_unbalanced_binomial()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)

  result <- der_diagnose(
    draws_all,
    y            = fix$y,
    X            = fix$X,
    group        = fix$group,
    weights      = fix$w,
    psu          = fix$psu,
    family       = "binomial",
    sigma_theta  = fix$sigma_theta_hat,
    beta_prior_sd = fix$beta_prior_sd,
    param_types  = fix$param_types,
    tau          = 1.2
  )

  out <- withr::with_output_sink(tempfile(), summary(result))
  expected_cols <- c("param_name", "param_type", "der",
                     "tier", "tier_label", "flagged", "action")
  expect_true(all(expected_cols %in% names(out)))
})

test_that("summary works on unclassified svyder object", {
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

  out <- withr::with_output_sink(tempfile(), summary(result))
  expect_true(is.data.frame(out))
  expect_true("param_name" %in% names(out))
  expect_true("der" %in% names(out))
})

test_that("summary row count matches parameter count", {
  fix <- make_unbalanced_binomial()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)

  result <- der_diagnose(
    draws_all,
    y            = fix$y,
    X            = fix$X,
    group        = fix$group,
    weights      = fix$w,
    psu          = fix$psu,
    family       = "binomial",
    sigma_theta  = fix$sigma_theta_hat,
    beta_prior_sd = fix$beta_prior_sd,
    param_types  = fix$param_types,
    tau          = 1.2
  )

  out <- withr::with_output_sink(tempfile(), summary(result))
  expect_equal(nrow(out), fix$d)
})
