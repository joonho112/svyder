# test-classes.R
# Tests for the S3 class infrastructure (new_svyder, validate_svyder, is.svyder)


# --- Helper to build a minimal valid svyder object for testing ---
.make_test_svyder <- function(d = 3, J = 2, M = 100) {

  p <- d - J

  der            <- setNames(runif(d, 0.5, 1.5), paste0("param", seq_len(d)))
  params         <- paste0("param", seq_len(d))
  H_obs          <- diag(d)
  J_cluster      <- diag(d)
  V_sand         <- diag(d)
  sigma_mcmc     <- diag(d)

  deff_j         <- rep(1.0, J)
  B_j            <- rep(0.8, J)

  classification <- data.frame(
    param_name = params,
    param_type = c(rep("fe_between", p), rep("re", J)),
    der        = as.numeric(der),
    tier       = c(rep("I-b", p), rep("II", J)),
    flagged    = rep(FALSE, d),
    action     = rep("retain", d),
    stringsAsFactors = FALSE
  )

  corrected_draws <- matrix(rnorm(M * d), nrow = M, ncol = d)
  original_draws  <- matrix(rnorm(M * d), nrow = M, ncol = d)
  scale_factors   <- rep(1.0, d)

  new_svyder(
    der              = der,
    params           = params,
    H_obs            = H_obs,
    J_cluster        = J_cluster,
    V_sand           = V_sand,
    sigma_mcmc       = sigma_mcmc,
    deff_j           = deff_j,
    B_j              = B_j,
    classification   = classification,
    tau              = 1.2,
    corrected_draws  = corrected_draws,
    scale_factors    = scale_factors,
    original_draws   = original_draws,
    call             = match.call(),
    family           = "gaussian",
    n_obs            = 100L,
    n_groups         = J,
    compute_time     = 0.5
  )
}


# ============================================================================
# new_svyder tests
# ============================================================================

test_that("new_svyder creates object with class 'svyder'", {
  obj <- .make_test_svyder()
  expect_s3_class(obj, "svyder")
})

test_that("new_svyder returns a list", {
  obj <- .make_test_svyder()
  expect_type(obj, "list")
})

test_that("new_svyder stores all fields", {
  obj <- .make_test_svyder()

  expected_fields <- c(
    "der", "params", "H_obs", "J_cluster", "V_sand", "sigma_mcmc",
    "deff_j", "B_j", "classification", "tau", "corrected_draws",
    "scale_factors", "original_draws", "call", "family",
    "n_obs", "n_groups", "compute_time"
  )

  for (f in expected_fields) {
    expect_true(f %in% names(obj),
                info = paste("Field", f, "should be present"))
  }
})

test_that("new_svyder preserves input values", {
  obj <- .make_test_svyder(d = 4, J = 2, M = 50)

  expect_equal(length(obj$der), 4)
  expect_equal(length(obj$params), 4)
  expect_equal(nrow(obj$H_obs), 4)
  expect_equal(ncol(obj$H_obs), 4)
  expect_equal(obj$tau, 1.2)
  expect_equal(obj$family, "gaussian")
  expect_equal(obj$n_obs, 100L)
  expect_equal(obj$n_groups, 2)
  expect_equal(nrow(obj$corrected_draws), 50)
})


# ============================================================================
# is.svyder tests
# ============================================================================

test_that("is.svyder returns TRUE for svyder objects", {
  obj <- .make_test_svyder()
  expect_true(is.svyder(obj))
})

test_that("is.svyder returns FALSE for non-svyder objects", {
  expect_false(is.svyder(list(a = 1)))
  expect_false(is.svyder(1:10))
  expect_false(is.svyder("hello"))
  expect_false(is.svyder(NULL))
  expect_false(is.svyder(data.frame(x = 1)))
})

test_that("is.svyder returns FALSE for unclassed list with same structure", {
  obj <- .make_test_svyder()
  unclass_obj <- unclass(obj)
  expect_false(is.svyder(unclass_obj))
})


# ============================================================================
# validate_svyder tests
# ============================================================================

test_that("validate_svyder passes for valid object", {
  obj <- .make_test_svyder()
  expect_invisible(validate_svyder(obj))
  expect_s3_class(validate_svyder(obj), "svyder")
})

test_that("validate_svyder catches non-svyder input", {
  expect_error(validate_svyder(list(a = 1)),
               "not of class 'svyder'")
})

test_that("validate_svyder catches mismatched DER length", {
  obj <- .make_test_svyder(d = 3, J = 2)
  # Make params length mismatch
  obj$params <- c("a", "b")
  expect_error(validate_svyder(obj),
               "does not match length of 'der'")
})

test_that("validate_svyder catches mismatched H_obs dimensions", {
  obj <- .make_test_svyder(d = 3, J = 2)
  obj$H_obs <- diag(4)
  expect_error(validate_svyder(obj),
               "must be a 3 x 3 matrix")
})

test_that("validate_svyder catches mismatched classification rows", {
  obj <- .make_test_svyder(d = 3, J = 2)
  obj$classification <- obj$classification[1:2, ]
  expect_error(validate_svyder(obj),
               "does not match length of 'der'")
})

test_that("validate_svyder catches mismatched scale_factors length", {
  obj <- .make_test_svyder(d = 3, J = 2)
  obj$scale_factors <- c(1.0, 1.0)
  expect_error(validate_svyder(obj),
               "does not match length of 'der'")
})

test_that("validate_svyder catches mismatched deff_j length", {
  obj <- .make_test_svyder(d = 3, J = 2)
  obj$deff_j <- 1.0  # length 1, should be 2

  expect_error(validate_svyder(obj),
               "does not match 'n_groups'")
})

test_that("validate_svyder catches mismatched B_j length", {
  obj <- .make_test_svyder(d = 3, J = 2)
  obj$B_j <- rep(0.8, 3)  # length 3, should be 2
  expect_error(validate_svyder(obj),
               "does not match 'n_groups'")
})

test_that("validate_svyder catches non-finite tau", {
  obj <- .make_test_svyder()
  obj$tau <- Inf
  expect_error(validate_svyder(obj),
               "finite numeric scalar")
})

test_that("validate_svyder catches missing required fields", {
  obj <- .make_test_svyder()
  obj$V_sand <- NULL
  expect_error(validate_svyder(obj),
               "Missing required fields")
})

test_that("validate_svyder catches mismatched corrected_draws columns", {
  obj <- .make_test_svyder(d = 3, J = 2, M = 50)
  obj$corrected_draws <- matrix(1, nrow = 50, ncol = 4)
  expect_error(validate_svyder(obj),
               "must be a matrix with 3 columns")
})
