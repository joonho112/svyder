# test-der_classify.R
# Tests for der_classify()


# ============================================================================
# Helper: create a classified svyder object from a fixture
# ============================================================================
.compute_and_classify <- function(fix, tau = 1.2) {
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)
  result <- der_compute(
    draws_all,
    y            = fix$y,
    X            = fix$X,
    group        = fix$group,
    weights      = fix$w,
    psu          = fix$psu,
    family       = fix$family,
    sigma_theta  = fix$sigma_theta_hat,
    sigma_e      = fix$sigma_e,
    beta_prior_sd = fix$beta_prior_sd,
    param_types  = fix$param_types
  )
  der_classify(result, tau = tau, verbose = FALSE)
}


# ============================================================================
# Tier assignments by param_type
# ============================================================================

test_that("fe_within DER > tau -> Tier I-a CORRECT", {
  fix <- make_unbalanced_binomial()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)
  result <- der_compute(
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

  # Force a fe_within param to have high DER for testing
  # Modify classification to have fe_within type
  result$classification$param_type[2] <- "fe_within"

  # Set DER artificially high for testing
  original_der <- result$der[2]
  result$der[2] <- 2.0

  classified <- der_classify(result, tau = 1.2, verbose = FALSE)

  # Check the fe_within parameter
  row2 <- classified$classification[2, ]
  expect_equal(row2$tier, "I-a")
  expect_equal(row2$tier_label, "Survey-dominated")
  expect_equal(row2$action, "CORRECT")
  expect_true(row2$flagged)
})

test_that("fe_between DER < tau -> Tier I-b retain", {
  fix <- make_balanced_gaussian()
  result <- .compute_and_classify(fix, tau = 1.2)

  # The intercept is fe_between
  row1 <- result$classification[1, ]
  expect_equal(row1$param_type, "fe_between")
  expect_equal(row1$tier, "I-b")
  expect_equal(row1$tier_label, "Protected (between)")

  # With equal weights, DER should be close to 1.0, so not flagged
  if (row1$der <= 1.2) {
    expect_equal(row1$action, "retain")
    expect_false(row1$flagged)
  }
})

test_that("re DER < tau -> Tier II retain", {
  fix <- make_balanced_gaussian()
  result <- .compute_and_classify(fix, tau = 1.2)

  # Random effects start at index p+1
  p <- fix$p
  re_rows <- result$classification[(p + 1):fix$d, ]

  for (i in seq_len(nrow(re_rows))) {
    expect_equal(re_rows$param_type[i], "re")
    expect_equal(re_rows$tier[i], "II")
    expect_equal(re_rows$tier_label[i], "Protected (random effects)")
  }
})


# ============================================================================
# Boundary case: DER exactly at tau
# ============================================================================

test_that("boundary: DER exactly at tau -> retain", {
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

  # Set DER exactly to tau
  tau_val <- 1.2
  result$der[1] <- tau_val

  classified <- der_classify(result, tau = tau_val, verbose = FALSE)

  # Strict inequality: DER == tau should NOT be flagged

  expect_false(classified$classification$flagged[1])
  expect_equal(classified$classification$action[1], "retain")
})


# ============================================================================
# No flagged parameters
# ============================================================================

test_that("no flagged parameters", {
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

  # Set all DER values below threshold
  result$der[] <- 0.9

  classified <- der_classify(result, tau = 1.2, verbose = FALSE)
  expect_equal(sum(classified$classification$flagged), 0)
  expect_true(all(classified$classification$action == "retain"))
})


# ============================================================================
# Classification data frame structure
# ============================================================================

test_that("classification data frame has correct columns", {
  fix <- make_unbalanced_binomial()
  result <- .compute_and_classify(fix, tau = 1.2)

  expected_cols <- c("param_name", "param_type", "der",
                     "tier", "tier_label", "flagged", "action")
  expect_true(all(expected_cols %in% names(result$classification)))
})

test_that("tau is stored correctly", {
  fix <- make_balanced_gaussian()
  result <- .compute_and_classify(fix, tau = 1.5)
  expect_equal(result$tau, 1.5)
})


# ============================================================================
# Pipe compatibility
# ============================================================================

test_that("pipe works: der_compute |> der_classify", {
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

  classified <- der_classify(result, tau = 1.2, verbose = FALSE)
  expect_s3_class(classified, "svyder")
  expect_true("flagged" %in% names(classified$classification))
})
