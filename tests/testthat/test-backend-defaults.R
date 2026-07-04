# test-backend-defaults.R
# Lightweight checks for backend method API defaults.

test_that("all der_compute methods default to diffuse beta prior", {
  expect_identical(formals(svyder:::der_compute.matrix)$beta_prior_sd, Inf)
  expect_identical(formals(svyder:::der_compute.brmsfit)$beta_prior_sd, Inf)
  expect_identical(formals(svyder:::der_compute.CmdStanMCMC)$beta_prior_sd, Inf)
  expect_identical(formals(svyder:::der_compute.stanreg)$beta_prior_sd, Inf)
})

test_that("matrix path records omitted and explicit beta_prior_sd targets", {
  fix <- make_balanced_gaussian()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)

  default_result <- der_compute(
    draws_all,
    y = fix$y, X = fix$X, group = fix$group,
    weights = fix$w, cluster = fix$psu,
    family = fix$family,
    sigma_theta = fix$sigma_theta_hat,
    sigma_e = fix$sigma_e,
    param_types = fix$param_types
  )

  finite_result <- der_compute(
    draws_all,
    y = fix$y, X = fix$X, group = fix$group,
    weights = fix$w, cluster = fix$psu,
    family = fix$family,
    sigma_theta = fix$sigma_theta_hat,
    sigma_e = fix$sigma_e,
    beta_prior_sd = 5,
    param_types = fix$param_types
  )

  expect_identical(default_result$target$beta_prior_sd, Inf)
  expect_identical(finite_result$target$beta_prior_sd, 5)
})
