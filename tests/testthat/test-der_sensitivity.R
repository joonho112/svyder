# test-der_sensitivity.R
# Tests for DER sensitivity analysis across threshold values

# --- Helper: build a svyder object from fixture ---
.fixture_to_svyder_sens <- function(fix) {
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


test_that("sensitivity returns correct columns", {
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder_sens(fix)
  sen <- der_sensitivity(sv)

  expect_s3_class(sen, "data.frame")
  expected_cols <- c("tau", "n_flagged", "pct_flagged", "flagged_params")
  expect_true(all(expected_cols %in% names(sen)))
})


test_that("n_flagged is monotone non-increasing", {
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder_sens(fix)
  sen <- der_sensitivity(sv, tau_range = seq(0.1, 5.0, by = 0.1))

  # n_flagged should be non-increasing as tau increases
  diffs <- diff(sen$n_flagged)
  expect_true(all(diffs <= 0))
})


test_that("tau = 0 flags everything", {
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder_sens(fix)

  # All DER values > 0 should be flagged at tau = 0
  # (DER > 0 is true for any non-degenerate posterior)
  sen <- der_sensitivity(sv, tau_range = 0)
  expect_equal(sen$n_flagged, fix$d)
})


test_that("tau = Inf flags nothing", {
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder_sens(fix)

  sen <- der_sensitivity(sv, tau_range = Inf)
  expect_equal(sen$n_flagged, 0L)
})


test_that("default tau range works", {
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder_sens(fix)
  sen <- der_sensitivity(sv)

  expect_equal(nrow(sen), length(seq(0.8, 2.0, by = 0.1)))
  expect_true(all(sen$pct_flagged >= 0 & sen$pct_flagged <= 1))
})


test_that("flagged_params is a list-column of character vectors", {
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder_sens(fix)
  sen <- der_sensitivity(sv)

  expect_true(is.list(sen$flagged_params))
  for (k in seq_len(nrow(sen))) {
    expect_true(is.character(sen$flagged_params[[k]]))
  }
})


test_that("flagged_params length matches n_flagged", {
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder_sens(fix)
  sen <- der_sensitivity(sv)

  for (k in seq_len(nrow(sen))) {
    expect_equal(length(sen$flagged_params[[k]]), sen$n_flagged[k])
  }
})


test_that("sensitivity rejects non-svyder input", {
  expect_error(der_sensitivity("not_svyder"))
  expect_error(der_sensitivity(list(a = 1)))
})


test_that("custom tau_range works", {
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder_sens(fix)
  custom_tau <- c(0.5, 1.0, 1.5, 2.0, 3.0)
  sen <- der_sensitivity(sv, tau_range = custom_tau)

  expect_equal(nrow(sen), length(custom_tau))
  expect_equal(sen$tau, sort(custom_tau))
})
