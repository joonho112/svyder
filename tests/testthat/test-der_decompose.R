# test-der_decompose.R
# Tests for the structural DER decomposition into constituent factors

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
    cluster = fix$psu,
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
  expected_cols <- c("param", "param_type", "der", "deff_used",
                     "B_used", "R_k", "kappa", "der_predicted")
  expect_true(all(expected_cols %in% names(dec)))
  expect_equal(nrow(dec), fix$d)
})


test_that("der_decompose values non-negative", {
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder(fix)
  dec <- der_decompose(sv)

  expect_true(all(dec$der >= 0))
  expect_true(all(dec$deff_used >= 0))
  expect_true(all(dec$B_used[!is.na(dec$B_used)] >= 0))
  expect_true(all(dec$der_predicted[!is.na(dec$der_predicted)] >= 0))
})


test_that("B_used in [0, 1] for RE rows and NA for FE rows", {
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder(fix)
  dec <- der_decompose(sv)

  re_rows <- dec$param_type == "re"
  fe_rows <- dec$param_type %in% c("fe_within", "fe_between")

  expect_true(all(dec$B_used[re_rows] >= 0 & dec$B_used[re_rows] <= 1))
  expect_true(all(is.na(dec$B_used[fe_rows])))
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


test_that("R_k is structural: intercept R_k equals information-weighted mean B", {
  # For an intercept column, the structural formula reduces to
  # R_1 = sum_j(B_j * a_j) / sum_j(a_j) with a_j = sum of working weights.
  fix <- make_balanced_gaussian()   # intercept-only, balanced, equal weights
  sv  <- .fixture_to_svyder(fix)
  dec <- der_decompose(sv)

  a_j <- as.numeric(tapply(sv$data$v, sv$data$group, sum))
  expected_R1 <- sum(sv$B_j * a_j) / sum(a_j)

  expect_equal(dec$R_k[1], expected_R1, tolerance = 1e-12)
  # For the balanced fixture this equals the analytical B = 5/6
  expect_equal(dec$R_k[1], fix$B_analytical, tolerance = 1e-10)
})


test_that("FE prediction is DEFF * (1 - R_k) with the global Kish DEFF", {
  fix <- make_unbalanced_binomial()
  sv  <- .fixture_to_svyder(fix)
  dec <- der_decompose(sv)

  w <- sv$data$weights
  deff_global <- length(w) * sum(w^2) / sum(w)^2

  fe_rows <- which(dec$param_type %in% c("fe_within", "fe_between"))
  expect_equal(dec$deff_used[fe_rows], rep(deff_global, length(fe_rows)),
               tolerance = 1e-12)
  expect_equal(dec$der_predicted[fe_rows],
               deff_global * (1 - dec$R_k[fe_rows]),
               tolerance = 1e-12)
})


test_that("RE prediction matches per-group B_j * DEFF_j * kappa_j", {
  fix <- make_unbalanced_binomial()
  sv  <- .fixture_to_svyder(fix)
  dec <- der_decompose(sv)

  re_rows <- which(dec$param_type == "re")
  kappa_j <- svyder:::.compute_kappa_j(sv$B_j, sv$n_groups)
  expected <- as.numeric(sv$B_j * sv$deff_j * kappa_j)

  expect_equal(dec$der_predicted[re_rows], expected, tolerance = 1e-12)
  expect_equal(dec$B_used[re_rows], as.numeric(sv$B_j), tolerance = 1e-12)
  expect_equal(dec$deff_used[re_rows], as.numeric(sv$deff_j), tolerance = 1e-12)
})


test_that("der_predicted is NOT back-solved from the observed DER", {
  # The structural prediction is a genuine approximation: with MCMC noise
  # it should not coincide with the observed DER to machine precision.
  fix <- make_unbalanced_binomial()
  sv  <- .fixture_to_svyder(fix)
  dec <- der_decompose(sv)

  fe_rows <- dec$param_type %in% c("fe_within", "fe_between")
  expect_false(isTRUE(all.equal(dec$der_predicted[fe_rows],
                                dec$der[fe_rows],
                                tolerance = 1e-10)))
})


test_that("structural prediction tracks the observed DER (mechanism check)", {
  # Balanced gaussian, equal weights: theory predicts intercept DER
  # ~ DEFF * (1 - B) = 1 * (1 - 5/6) = 1/6. The observed DER should be in
  # the same regime (loose band; MCMC and realized-data noise).
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder(fix)
  dec <- der_decompose(sv)

  pred <- dec$der_predicted[1]
  obs  <- dec$der[1]
  expect_equal(pred, 1 - fix$B_analytical, tolerance = 1e-8)
  expect_gt(obs, pred / 3)
  expect_lt(obs, pred * 3)
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


test_that("der_decompose errors on objects without stored data slots", {
  fix <- make_balanced_gaussian()
  sv  <- .fixture_to_svyder(fix)
  sv$data <- NULL   # mimic an object created by svyder < 0.2.0
  expect_error(der_decompose(sv), "no stored data")
})


test_that("der_decompose rejects non-svyder input", {
  expect_error(der_decompose("not_svyder"))
  expect_error(der_decompose(list(a = 1)))
})
