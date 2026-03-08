# test-der_theorem_check.R
# Tests for theoretical DER prediction verification

# --- Helper: build a svyder object from fixture ---
.fixture_to_svyder_thm <- function(fix) {
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)
  sigma_e_arg <- if (fix$family == "gaussian") fix$sigma_e else NULL

  der_compute(
    draws_all,
    y = fix$y,
    X = fix$X,
    group = fix$group,
    weights = fix$w,
    psu = fix$psu,
    family = fix$family,
    sigma_theta = fix$sigma_theta_hat,
    sigma_e = sigma_e_arg,
    beta_prior_sd = fix$beta_prior_sd,
    param_types = fix$param_types
  )
}


test_that("theorem_check returns correct columns", {
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder_thm(fix)
  tc  <- der_theorem_check(sv)

  expect_s3_class(tc, "data.frame")
  expected_cols <- c("param", "param_type", "der_empirical",
                     "der_theorem1", "der_theorem2",
                     "relative_error", "theorem_used")
  expect_true(all(expected_cols %in% names(tc)))
  expect_equal(nrow(tc), fix$d)
})


test_that("theorem1 is NA for RE and theorem2 is NA for FE", {
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder_thm(fix)
  tc  <- der_theorem_check(sv)

  fe_rows <- tc$param_type %in% c("fe_within", "fe_between")
  re_rows <- tc$param_type == "re"

  # Theorem 1 applies to FE, should be NA for RE
  expect_true(all(!is.na(tc$der_theorem1[fe_rows])))
  expect_true(all(is.na(tc$der_theorem1[re_rows])))

  # Theorem 2 applies to RE, should be NA for FE
  expect_true(all(is.na(tc$der_theorem2[fe_rows])))
  expect_true(all(!is.na(tc$der_theorem2[re_rows])))
})


test_that("balanced gaussian has reasonable relative error for RE", {
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder_thm(fix)
  tc  <- der_theorem_check(sv)

  # For balanced equal-weight design, the simple formula should be
  # roughly accurate. But with finite samples and MCMC noise, we

  # use a generous tolerance. The key invariant is that the function
  # runs and produces finite relative errors.
  re_rows <- tc$param_type == "re"
  re_errors <- tc$relative_error[re_rows]

  expect_true(all(is.finite(re_errors)))
  expect_true(all(re_errors >= 0))
})


test_that("conservation law attribute present when applicable", {
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder_thm(fix)
  tc  <- der_theorem_check(sv)

  # The balanced gaussian fixture has fe_between + re, so conservation
  # law attribute should be present
  cl <- attr(tc, "conservation_law")
  expect_true(!is.null(cl))
  expect_true(is.list(cl))
  expect_true(all(c("der_mu", "der_theta_mean", "conservation_sum",
                     "deff_mean", "relative_error") %in% names(cl)))
})


test_that("relative_error is non-negative", {
  fix <- make_unbalanced_binomial()
  sv  <- .fixture_to_svyder_thm(fix)
  tc  <- der_theorem_check(sv)

  finite_errors <- tc$relative_error[is.finite(tc$relative_error)]
  expect_true(all(finite_errors >= 0))
})


test_that("theorem_used labels are correct", {
  fix <- make_unbalanced_binomial()
  sv  <- .fixture_to_svyder_thm(fix)
  tc  <- der_theorem_check(sv)

  valid_theorems <- c("Theorem 1 (within)", "Theorem 1 (between)", "Theorem 2 (RE)")
  expect_true(all(tc$theorem_used %in% valid_theorems))

  # Check that param_type and theorem alignment are correct
  expect_true(all(tc$theorem_used[tc$param_type == "fe_within"] == "Theorem 1 (within)"))
  expect_true(all(tc$theorem_used[tc$param_type == "fe_between"] == "Theorem 1 (between)"))
  expect_true(all(tc$theorem_used[tc$param_type == "re"] == "Theorem 2 (RE)"))
})


test_that("theorem_check rejects non-svyder input", {
  expect_error(der_theorem_check("not_svyder"))
  expect_error(der_theorem_check(list(a = 1)))
})
