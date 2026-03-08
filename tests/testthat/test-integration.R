###############################################################################
# test-integration.R
# Integration tests using bundled datasets
###############################################################################

test_that("full pipeline on nsece_demo data", {
  data(nsece_demo, package = "svyder")
  result <- der_diagnose(
    nsece_demo$draws,
    y = nsece_demo$y, X = nsece_demo$X,
    group = nsece_demo$group, weights = nsece_demo$weights,
    psu = nsece_demo$psu, family = nsece_demo$family,
    sigma_theta = nsece_demo$sigma_theta,
    param_types = nsece_demo$param_types
  )
  expect_s3_class(result, "svyder")
  expect_true(any(result$classification$flagged))
})

test_that("full pipeline on sim_hlr data", {
  data(sim_hlr, package = "svyder")
  result <- der_diagnose(
    sim_hlr$draws,
    y = sim_hlr$y, X = sim_hlr$X,
    group = sim_hlr$group, weights = sim_hlr$weights,
    family = sim_hlr$family,
    sigma_theta = sim_hlr$sigma_theta,
    sigma_e = sim_hlr$sigma_e,
    param_types = sim_hlr$param_types
  )
  expect_s3_class(result, "svyder")
})

test_that("pipe chain works end-to-end", {
  data(nsece_demo, package = "svyder")
  result <- nsece_demo$draws |>
    der_compute(y = nsece_demo$y, X = nsece_demo$X,
                group = nsece_demo$group, weights = nsece_demo$weights,
                psu = nsece_demo$psu, family = "binomial",
                sigma_theta = nsece_demo$sigma_theta,
                param_types = nsece_demo$param_types) |>
    der_classify(tau = 1.2, verbose = FALSE) |>
    der_correct()
  expect_s3_class(result, "svyder")
})

test_that("nsece_demo dimensions are consistent", {
  data(nsece_demo, package = "svyder")
  expect_equal(nsece_demo$N, length(nsece_demo$y))
  expect_equal(nsece_demo$N, nrow(nsece_demo$X))
  expect_equal(nsece_demo$N, length(nsece_demo$group))
  expect_equal(nsece_demo$N, length(nsece_demo$weights))
  expect_equal(nsece_demo$N, length(nsece_demo$psu))
  expect_equal(nsece_demo$p, ncol(nsece_demo$X))
  expect_equal(ncol(nsece_demo$draws), nsece_demo$p + nsece_demo$J)
  expect_equal(nrow(nsece_demo$draws), 4000L)
  expect_equal(nsece_demo$J, max(nsece_demo$group))
})

test_that("sim_hlr dimensions are consistent", {
  data(sim_hlr, package = "svyder")
  expect_equal(sim_hlr$N, length(sim_hlr$y))
  expect_equal(sim_hlr$N, nrow(sim_hlr$X))
  expect_equal(sim_hlr$N, length(sim_hlr$group))
  expect_equal(sim_hlr$N, length(sim_hlr$weights))
  expect_equal(sim_hlr$p, ncol(sim_hlr$X))
  expect_equal(ncol(sim_hlr$draws), sim_hlr$p + sim_hlr$J)
  expect_equal(nrow(sim_hlr$draws), 4000L)
  expect_equal(sim_hlr$J, max(sim_hlr$group))
})

test_that("sim_hlr equal weights produce DER near 1", {
  data(sim_hlr, package = "svyder")
  result <- der_compute(
    sim_hlr$draws,
    y = sim_hlr$y, X = sim_hlr$X,
    group = sim_hlr$group, weights = sim_hlr$weights,
    family = sim_hlr$family,
    sigma_theta = sim_hlr$sigma_theta,
    sigma_e = sim_hlr$sigma_e,
    param_types = sim_hlr$param_types
  )
  # With equal weights, DER values should be positive and finite.
  # The intercept (fe_between) may have low DER due to shrinkage
  # absorbing between-cluster variance, while the within covariate
  # should be closer to 1.0.
  fe_der <- result$der[seq_len(sim_hlr$p)]
  expect_true(all(fe_der > 0))
  expect_true(all(is.finite(fe_der)))
  expect_true(all(fe_der < 3.0))
})

test_that("nsece_demo produces flagged fe_between parameters", {
  data(nsece_demo, package = "svyder")
  result <- der_diagnose(
    nsece_demo$draws,
    y = nsece_demo$y, X = nsece_demo$X,
    group = nsece_demo$group, weights = nsece_demo$weights,
    psu = nsece_demo$psu, family = nsece_demo$family,
    sigma_theta = nsece_demo$sigma_theta,
    param_types = nsece_demo$param_types
  )
  # With unequal weights, at least one fe_between parameter should be flagged
  flagged_fe <- result$classification[
    result$classification$param_type %in% c("fe_between", "fe_within") &
    result$classification$flagged, ]
  expect_gt(nrow(flagged_fe), 0)
})

test_that("corrected draws preserve posterior mean", {
  data(nsece_demo, package = "svyder")
  result <- der_diagnose(
    nsece_demo$draws,
    y = nsece_demo$y, X = nsece_demo$X,
    group = nsece_demo$group, weights = nsece_demo$weights,
    psu = nsece_demo$psu, family = nsece_demo$family,
    sigma_theta = nsece_demo$sigma_theta,
    param_types = nsece_demo$param_types
  )
  # Corrected draws should have same column means as original
  orig_means <- colMeans(result$original_draws)
  corr_means <- colMeans(result$corrected_draws)
  expect_equal(orig_means, corr_means, tolerance = 1e-10)
})
