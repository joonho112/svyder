# test-der_compute.R
# Tests for the core der_compute() function


# ============================================================================
# Helper: run der_compute on a fixture list
# ============================================================================
.run_der_compute <- function(fix) {
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)
  der_compute(
    draws_all,
    y            = fix$y,
    X            = fix$X,
    group        = fix$group,
    weights      = fix$w,
    cluster      = fix$psu,
    family       = fix$family,
    sigma_theta  = fix$sigma_theta_hat,
    sigma_e      = fix$sigma_e,
    beta_prior_sd = fix$beta_prior_sd,
    param_types  = fix$param_types

  )
}


# ============================================================================
# Basic functionality
# ============================================================================

test_that("der_compute returns svyder class", {
  fix <- make_balanced_gaussian()
  result <- .run_der_compute(fix)
  expect_s3_class(result, "svyder")
})

test_that("DER values are positive", {
  fix <- make_balanced_gaussian()
  result <- .run_der_compute(fix)
  expect_true(all(result$der > 0))
})

test_that("DER length matches parameters", {
  fix <- make_balanced_gaussian()
  result <- .run_der_compute(fix)
  expected_d <- fix$p + fix$J
  expect_equal(length(result$der), expected_d)
  expect_equal(length(result$params), expected_d)
})

test_that("pipe-friendly: returns svyder", {
  fix <- make_balanced_gaussian()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)
  result <- der_compute(
    draws_all,
    y = fix$y, X = fix$X, group = fix$group,
    weights = fix$w, cluster = fix$psu,
    family = "gaussian", sigma_theta = fix$sigma_theta_hat,
    sigma_e = fix$sigma_e, beta_prior_sd = fix$beta_prior_sd,
    param_types = fix$param_types
  )
  expect_s3_class(result, "svyder")
  # Can pipe into der_classify
  result2 <- der_classify(result, tau = 1.2, verbose = FALSE)
  expect_s3_class(result2, "svyder")
})

test_that("errors on unknown class", {
  expect_error(
    der_compute("not_a_matrix"),
    "does not know how to handle class"
  )
  expect_error(
    der_compute(data.frame(x = 1:10)),
    "does not know how to handle class"
  )
})


# ============================================================================
# Regression test: matches standalone compute_der()
# ============================================================================

test_that("matches standalone compute_der()", {
  standalone_path <- file.path(
    "/Users/joonholee/Documents/2026-03-03_Design Effect Ratios Paper",
    "dev/scripts/standalone/compute_der.R"
  )
  skip_if_not(file.exists(standalone_path),
              "Standalone compute_der.R not found")

  source(standalone_path, local = TRUE)

  # Use binomial fixture (standalone only supports binomial)
  fix <- make_unbalanced_binomial()

  # --- Run standalone ---
  standalone_result <- compute_der(
    beta_hat         = fix$beta_hat,
    theta_hat        = fix$theta_hat,
    sigma_theta_hat  = fix$sigma_theta_hat,
    draws_beta       = fix$draws_beta,
    draws_theta      = fix$draws_theta,
    y                = fix$y,
    X                = fix$X,
    group            = fix$group,
    w                = fix$w,
    psu              = fix$psu,
    beta_prior_sd    = fix$beta_prior_sd,
    tau              = 1.2,
    param_types      = fix$param_types
  )

  # --- Run svyder (legacy meat flags reproduce the v1 pipeline) ---
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)
  svyder_result <- der_compute(
    draws_all,
    y            = fix$y,
    X            = fix$X,
    group        = fix$group,
    weights      = fix$w,
    cluster      = fix$psu,
    family       = "binomial",
    sigma_theta  = fix$sigma_theta_hat,
    center_meat  = FALSE,
    df_correct   = FALSE,
    beta_prior_sd = fix$beta_prior_sd,
    param_types  = fix$param_types
  )

  # DER values should be identical (same math path)
  # Note: der_compute uses colMeans(draws) as point estimates while
  # standalone uses fixed beta_hat/theta_hat. We need to account for this.
  # The standalone uses the true values as point_est, while der_compute
  # uses colMeans of draws. Since draws are generated from the true values,
  # the colMeans will be close but not identical. So we compare with a
  # small tolerance rather than exact equality.
  expect_equal(
    as.numeric(svyder_result$der),
    as.numeric(standalone_result$der),
    tolerance = 0.05,
    info = "DER values should closely match standalone"
  )

  # Sandwich matrices should also be close
  expect_equal(
    diag(svyder_result$V_sand),
    diag(standalone_result$V_sand),
    tolerance = 0.05
  )
})


# ============================================================================
# Gaussian fixture
# ============================================================================

test_that("balanced gaussian works end-to-end", {
  fix <- make_balanced_gaussian()
  result <- .run_der_compute(fix)

  # Should have correct dimensions
  expect_equal(result$n_obs, fix$N)
  expect_equal(result$n_groups, fix$J)
  expect_equal(result$family, "gaussian")

  # DER for equal weights should be positive and finite
  # Some random effects may have small DER due to strong shrinkage
  expect_true(all(result$der > 0))
  expect_true(all(is.finite(result$der)))
})


# ============================================================================
# Binomial fixture
# ============================================================================

test_that("binomial fixture works end-to-end", {
  fix <- make_unbalanced_binomial()
  result <- .run_der_compute(fix)

  expect_s3_class(result, "svyder")
  expect_equal(result$family, "binomial")
  expect_equal(result$n_obs, fix$N)
  expect_equal(result$n_groups, fix$J)
  expect_equal(length(result$der), fix$d)
  expect_true(all(result$der > 0))
})


# ============================================================================
# Minimal J=2 fixture
# ============================================================================

test_that("minimal J=2 fixture works", {
  fix <- make_minimal_j2()
  result <- .run_der_compute(fix)

  expect_s3_class(result, "svyder")
  expect_equal(result$n_groups, 2L)
  expect_equal(length(result$deff_j), 2L)
  expect_equal(length(result$B_j), 2L)
})


# ============================================================================
# Stored components
# ============================================================================

test_that("all matrix components have correct dimensions", {
  fix <- make_unbalanced_binomial()
  result <- .run_der_compute(fix)
  d <- fix$d

  expect_equal(dim(result$H_obs), c(d, d))
  expect_equal(dim(result$J_cluster), c(d, d))
  expect_equal(dim(result$V_sand), c(d, d))
  expect_equal(dim(result$sigma_mcmc), c(d, d))
  expect_equal(dim(result$original_draws), c(fix$M, d))
})

test_that("sigma_mcmc is symmetric positive definite", {
  fix <- make_balanced_gaussian()
  result <- .run_der_compute(fix)

  # Symmetric
  expect_equal(result$sigma_mcmc, t(result$sigma_mcmc))

  # Positive definite (all eigenvalues positive)
  eigenvals <- eigen(result$sigma_mcmc, only.values = TRUE)$values
  expect_true(all(eigenvals > 0))
})

test_that("V_sand is symmetric", {
  fix <- make_balanced_gaussian()
  result <- .run_der_compute(fix)
  expect_equal(result$V_sand, t(result$V_sand))
})


# ============================================================================
# Declared variance target: cluster is required
# ============================================================================

test_that("missing cluster errors informatively", {
  fix <- make_balanced_gaussian()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)
  expect_error(
    der_compute(
      draws_all,
      y = fix$y, X = fix$X, group = fix$group,
      weights = fix$w,
      family = "gaussian", sigma_theta = fix$sigma_theta_hat,
      sigma_e = fix$sigma_e, param_types = fix$param_types
    ),
    "'cluster' must be supplied explicitly"
  )
})

test_that("psu is a deprecated alias for cluster", {
  fix <- make_balanced_gaussian()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)

  result_cluster <- .run_der_compute(fix)

  expect_warning(
    result_psu <- der_compute(
      draws_all,
      y = fix$y, X = fix$X, group = fix$group,
      weights = fix$w, psu = fix$psu,
      family = "gaussian", sigma_theta = fix$sigma_theta_hat,
      sigma_e = fix$sigma_e, beta_prior_sd = fix$beta_prior_sd,
      param_types = fix$param_types
    ),
    "'psu' is deprecated"
  )
  expect_equal(as.numeric(result_psu$der), as.numeric(result_cluster$der),
               tolerance = 1e-12)
})

test_that("supplying both cluster and psu errors", {
  fix <- make_balanced_gaussian()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)
  expect_error(
    der_compute(
      draws_all,
      y = fix$y, X = fix$X, group = fix$group,
      weights = fix$w, cluster = fix$psu, psu = fix$psu,
      family = "gaussian", sigma_theta = fix$sigma_theta_hat,
      sigma_e = fix$sigma_e, param_types = fix$param_types
    ),
    "not both"
  )
})


# ============================================================================
# param_types inference
# ============================================================================

test_that("param_types = NULL emits a message and infers types", {
  fix <- make_unbalanced_binomial()  # X = [intercept, x1 (within-varying)]
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)

  expect_message(
    result <- der_compute(
      draws_all,
      y = fix$y, X = fix$X, group = fix$group,
      weights = fix$w, cluster = fix$psu,
      family = "binomial", sigma_theta = fix$sigma_theta_hat,
      beta_prior_sd = fix$beta_prior_sd
    ),
    "param_types inferred"
  )

  inferred <- result$classification$param_type[seq_len(fix$p)]
  expect_equal(inferred[1], "fe_between")  # intercept: constant within groups
  expect_equal(inferred[2], "fe_within")   # x1: varies within groups
})


# ============================================================================
# Weight normalization convention
# ============================================================================

test_that("normalize = 'unit_mean' is recorded and weights rescaled to mean 1", {
  fix <- make_unbalanced_binomial()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)

  # Rescale the fixture weights so raw mean != 1
  w_raw <- fix$w * 7

  result <- der_compute(
    draws_all,
    y = fix$y, X = fix$X, group = fix$group,
    weights = w_raw, cluster = fix$psu,
    family = "binomial", sigma_theta = fix$sigma_theta_hat,
    normalize = "unit_mean",
    beta_prior_sd = fix$beta_prior_sd,
    param_types = fix$param_types
  )

  expect_equal(result$target$normalize, "unit_mean")
  expect_equal(result$target$weights_raw_mean, mean(w_raw), tolerance = 1e-12)
  expect_equal(mean(result$data$weights), 1.0, tolerance = 1e-12)

  # DER must be invariant to overall weight rescaling under unit_mean
  result_base <- der_compute(
    draws_all,
    y = fix$y, X = fix$X, group = fix$group,
    weights = fix$w, cluster = fix$psu,
    family = "binomial", sigma_theta = fix$sigma_theta_hat,
    normalize = "unit_mean",
    beta_prior_sd = fix$beta_prior_sd,
    param_types = fix$param_types
  )
  expect_equal(as.numeric(result$der), as.numeric(result_base$der),
               tolerance = 1e-10)
})

test_that("normalize = 'none' uses weights as supplied", {
  fix <- make_unbalanced_binomial()
  draws_all <- cbind(fix$draws_beta, fix$draws_theta)

  result <- der_compute(
    draws_all,
    y = fix$y, X = fix$X, group = fix$group,
    weights = fix$w, cluster = fix$psu,
    family = "binomial", sigma_theta = fix$sigma_theta_hat,
    normalize = "none",
    beta_prior_sd = fix$beta_prior_sd,
    param_types = fix$param_types
  )

  expect_equal(result$target$normalize, "none")
  expect_equal(result$data$weights, fix$w, tolerance = 1e-15)
})


# ============================================================================
# Target metadata and Tier III exclusion table
# ============================================================================

test_that("target metadata records the declared variance target", {
  fix <- make_balanced_gaussian()
  result <- .run_der_compute(fix)

  tg <- result$target
  expect_type(tg, "list")
  expect_equal(tg$cluster_n, length(unique(fix$psu)))
  expect_true(tg$cluster_is_group)  # fixture uses psu = group
  expect_equal(tg$strata_n, 1L)
  expect_true(tg$center)
  expect_true(tg$df_correct)
  expect_equal(tg$normalize, "unit_mean")
  expect_equal(tg$sigma_theta, fix$sigma_theta_hat)
})

test_that("Tier III exclusion table lists hyperparameters", {
  fix_g <- make_balanced_gaussian()
  result_g <- .run_der_compute(fix_g)
  expect_s3_class(result_g$excluded, "data.frame")
  expect_setequal(result_g$excluded$param, c("sigma_theta", "sigma_e"))
  expect_true(all(result_g$excluded$tier == "III"))
  reason_theta <- result_g$excluded$reason[result_g$excluded$param == "sigma_theta"]
  reason_sigma_e <- result_g$excluded$reason[result_g$excluded$param == "sigma_e"]
  expect_match(reason_theta, "prior hyperparameter", fixed = TRUE)
  expect_match(reason_theta, "data-level phi=(beta, theta) score block", fixed = TRUE)
  expect_match(reason_sigma_e, "gaussian dispersion plug-in", fixed = TRUE)
  expect_match(reason_sigma_e, "location-block sandwich", fixed = TRUE)
  expect_false(identical(reason_theta, reason_sigma_e))

  fix_b <- make_unbalanced_binomial()
  result_b <- .run_der_compute(fix_b)
  expect_equal(result_b$excluded$param, "sigma_theta")
  expect_match(result_b$excluded$reason, "prior hyperparameter", fixed = TRUE)
})

test_that("data slots are stored for re-targeting", {
  fix <- make_unbalanced_binomial()
  result <- .run_der_compute(fix)

  expect_type(result$data, "list")
  expect_true(all(c("X", "group", "weights", "cluster", "strata", "v", "r")
                  %in% names(result$data)))
  expect_equal(nrow(result$data$X), fix$N)
  expect_equal(length(result$data$r), fix$N)
})
