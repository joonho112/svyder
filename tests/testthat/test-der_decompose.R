# test-der_decompose.R
# Tests for DER decomposition into constituent factors

# --- Helper: build a svyder object from fixture ---
.fixture_to_svyder <- function(fix) {
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


test_that("der_decompose returns correct columns", {
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder(fix)
  dec <- der_decompose(sv)

  expect_s3_class(dec, "data.frame")
  expected_cols <- c("param", "param_type", "der", "deff_mean",
                     "B_mean", "R_k", "kappa", "der_predicted")
  expect_true(all(expected_cols %in% names(dec)))
  expect_equal(nrow(dec), fix$d)
})


test_that("der_decompose values non-negative", {
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder(fix)
  dec <- der_decompose(sv)

  expect_true(all(dec$der >= 0))
  expect_true(all(dec$deff_mean >= 0))
  expect_true(all(dec$B_mean >= 0))
  expect_true(all(dec$der_predicted[!is.na(dec$der_predicted)] >= 0))
})


test_that("B values in [0, 1]", {
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder(fix)
  dec <- der_decompose(sv)

  expect_true(all(dec$B_mean >= 0 & dec$B_mean <= 1))
})


test_that("R_k values in [0, 1] for FE", {
  fix <- make_unbalanced_binomial()
  sv  <- .fixture_to_svyder(fix)
  dec <- der_decompose(sv)

  fe_rows <- dec$param_type %in% c("fe_within", "fe_between")
  rk_vals <- dec$R_k[fe_rows]
  expect_true(all(!is.na(rk_vals)))
  expect_true(all(rk_vals >= 0 & rk_vals <= 1))
})


test_that("kappa values in [0, 1] for RE", {
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder(fix)
  dec <- der_decompose(sv)

  re_rows    <- dec$param_type == "re"
  kappa_vals <- dec$kappa[re_rows]
  expect_true(all(!is.na(kappa_vals)))
  expect_true(all(kappa_vals >= 0 & kappa_vals <= 1))
})


test_that("decomposition consistent with DER", {
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder(fix)
  dec <- der_decompose(sv)

  # For FE parameters: der_predicted = deff_mean * (1 - R_k)
  # should match der since R_k is back-solved
  fe_rows <- dec$param_type %in% c("fe_within", "fe_between")
  if (any(fe_rows)) {
    expect_equal(
      dec$der_predicted[fe_rows],
      dec$der[fe_rows],
      tolerance = 1e-10
    )
  }
})


test_that("R_k is NA for random effects", {
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder(fix)
  dec <- der_decompose(sv)

  re_rows <- dec$param_type == "re"
  expect_true(all(is.na(dec$R_k[re_rows])))
})


test_that("kappa is NA for fixed effects", {
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder(fix)
  dec <- der_decompose(sv)

  fe_rows <- dec$param_type %in% c("fe_within", "fe_between")
  expect_true(all(is.na(dec$kappa[fe_rows])))
})


test_that("der_decompose rejects non-svyder input", {
  expect_error(der_decompose("not_svyder"))
  expect_error(der_decompose(list(a = 1)))
})
