# test-invariance.R
# Units-invariance of the DER
# ---------------------------------------------------------------------------
# The gaussian score is w * (y - mu) / sigma_e^2, matching the bread's
# working weights w / sigma_e^2. Under a change of response units
# (y, sigma_e, sigma_theta, draws) -> (c*y, c*sigma_e, c*sigma_theta,
# c*draws), the sandwich and the MCMC covariance both scale by c^2, so
# the DER must be exactly invariant. The pre-0.2.0 score w * (y - mu)
# made the DER scale like sigma_e^4; this file guards the fix.


# Tiny gaussian fixture: J=5 groups, n_j=20, p=2 (intercept + within x),
# unequal weights, known sigma_e = 2.
.make_units_fixture <- function(seed = 202) {
  set.seed(seed)

  J   <- 5L
  n_j <- 20L
  N   <- J * n_j
  p   <- 2L
  d   <- p + J
  S   <- 500L

  sigma_e     <- 2.0
  sigma_theta <- 0.7

  group <- rep(seq_len(J), each = n_j)
  x1    <- rnorm(N)
  X     <- cbind(intercept = 1, x1 = x1)

  w <- exp(rnorm(N, sd = 0.4))

  beta_true  <- c(1.0, -0.5)
  theta_true <- rnorm(J, sd = sigma_theta)
  y <- as.numeric(X %*% beta_true + theta_true[group] +
                    rnorm(N, sd = sigma_e))

  # Fake posterior draws centered at the true values
  point <- c(beta_true, theta_true)
  draws <- sweep(matrix(rnorm(S * d, sd = 0.1), S, d), 2L, point, "+")

  list(
    y = y, X = X, group = group, w = w, psu = group,
    draws = draws, sigma_e = sigma_e, sigma_theta = sigma_theta,
    param_types = c("fe_between", "fe_within"),
    N = N, J = J, p = p, d = d
  )
}

.compute_units_der <- function(fx, scale = 1) {
  der_compute(
    fx$draws * scale,
    y           = fx$y * scale,
    X           = fx$X,
    group       = fx$group,
    weights     = fx$w,
    cluster     = fx$psu,
    family      = "gaussian",
    sigma_theta = fx$sigma_theta * scale,
    sigma_e     = fx$sigma_e * scale,
    param_types = fx$param_types
  )
}


test_that("gaussian DER is invariant to a change of response units", {
  fx <- .make_units_fixture()

  base   <- .compute_units_der(fx, scale = 1)
  scaled <- .compute_units_der(fx, scale = 3)

  expect_equal(as.numeric(scaled$der), as.numeric(base$der),
               tolerance = 1e-10)

  # The invariance holds for the whole sandwich, up to the c^2 factor
  expect_equal(scaled$V_sand, 9 * base$V_sand, tolerance = 1e-8)
  expect_equal(scaled$sigma_mcmc, 9 * base$sigma_mcmc, tolerance = 1e-8)
})


test_that("gaussian units-invariance holds under legacy meat flags too", {
  fx <- .make_units_fixture()

  base <- der_compute(
    fx$draws,
    y = fx$y, X = fx$X, group = fx$group, weights = fx$w,
    cluster = fx$psu, family = "gaussian",
    sigma_theta = fx$sigma_theta, sigma_e = fx$sigma_e,
    center_meat = FALSE, df_correct = FALSE,
    param_types = fx$param_types
  )
  scaled <- der_compute(
    fx$draws * 3,
    y = fx$y * 3, X = fx$X, group = fx$group, weights = fx$w,
    cluster = fx$psu, family = "gaussian",
    sigma_theta = fx$sigma_theta * 3, sigma_e = fx$sigma_e * 3,
    center_meat = FALSE, df_correct = FALSE,
    param_types = fx$param_types
  )

  expect_equal(as.numeric(scaled$der), as.numeric(base$der),
               tolerance = 1e-10)
})


test_that("binomial DER ignores sigma_e (no-op sanity)", {
  fix <- make_unbalanced_binomial()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)

  args <- list(
    draws_all,
    y = fix$y, X = fix$X, group = fix$group,
    weights = fix$w, cluster = fix$psu,
    family = "binomial", sigma_theta = fix$sigma_theta_hat,
    beta_prior_sd = fix$beta_prior_sd, param_types = fix$param_types
  )

  res_null    <- do.call(der_compute, args)
  res_sigma_e <- do.call(der_compute, c(args, list(sigma_e = 5)))

  expect_identical(as.numeric(res_null$der), as.numeric(res_sigma_e$der))
  expect_identical(res_null$V_sand, res_sigma_e$V_sand)
})
