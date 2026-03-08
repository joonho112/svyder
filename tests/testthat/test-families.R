###############################################################################
# test-families.R
# Tests for family-specific computations (R/families.R)
###############################################################################

# --- Binomial working weights ---

test_that("binomial working weights correct at mu=0.5", {
  # w * 0.5 * 0.5 = w * 0.25
  mu <- rep(0.5, 10)
  w <- rep(1, 10)
  expect_equal(.working_weights_binomial(mu, w), rep(0.25, 10))
})

test_that("binomial working weights correct at mu=0.1", {
  mu <- rep(0.1, 5)
  w <- rep(2, 5)
  expect_equal(.working_weights_binomial(mu, w), rep(2 * 0.1 * 0.9, 5))
})

test_that("binomial working weights correct at boundaries", {
  # At mu=0 and mu=1, working weights should be 0
  mu <- c(0, 1)
  w <- c(1, 1)
  expect_equal(.working_weights_binomial(mu, w), c(0, 0))
})

test_that("binomial working weights scale linearly with survey weights", {
  mu <- c(0.2, 0.5, 0.8)
  w1 <- rep(1, 3)
  w2 <- rep(3, 3)
  expect_equal(.working_weights_binomial(mu, w2),
               3 * .working_weights_binomial(mu, w1))
})


# --- Gaussian working weights ---

test_that("gaussian working weights correct", {
  w <- c(1, 2, 3)
  sigma_e <- 2.0
  expect_equal(.working_weights_gaussian(w, sigma_e), w / 4)
})

test_that("gaussian working weights with sigma_e=1", {
  w <- c(1, 2, 3)
  expect_equal(.working_weights_gaussian(w, 1.0), w)
})

test_that("gaussian working weights error on NULL sigma_e", {
  expect_error(.working_weights_gaussian(c(1, 2), NULL),
               "sigma_e")
})

test_that("gaussian working weights error on non-positive sigma_e", {
  expect_error(.working_weights_gaussian(c(1, 2), 0),
               "sigma_e")
  expect_error(.working_weights_gaussian(c(1, 2), -1),
               "sigma_e")
})


# --- Dispatcher ---

test_that("working_weights dispatcher routes to binomial correctly", {
  mu <- rep(0.5, 5)
  w <- rep(1, 5)
  expect_equal(.working_weights("binomial", mu, w),
               .working_weights_binomial(mu, w))
})

test_that("working_weights dispatcher routes to gaussian correctly", {
  w <- c(1, 2, 3)
  sigma_e <- 1.5
  expect_equal(.working_weights("gaussian", rep(0.5, 3), w, sigma_e = sigma_e),
               .working_weights_gaussian(w, sigma_e))
})

test_that("working_weights errors on unknown family", {
  expect_error(.working_weights("poisson", rep(0.5, 5), rep(1, 5)),
               "Unknown family")
})


# --- compute_mu ---

test_that("compute_mu binomial is logistic", {
  eta <- c(-2, 0, 2)
  mu <- .compute_mu("binomial", eta)
  expect_equal(mu, 1 / (1 + exp(-eta)))
})

test_that("compute_mu binomial at eta=0 gives 0.5", {
  expect_equal(.compute_mu("binomial", 0), 0.5)
})

test_that("compute_mu gaussian is identity", {
  eta <- c(-2, 0, 2)
  expect_equal(.compute_mu("gaussian", eta), eta)
})

test_that("compute_mu errors on unknown family", {
  expect_error(.compute_mu("poisson", c(0, 1)),
               "Unknown family")
})


# --- compute_residuals ---

test_that("residuals are w * (y - mu)", {
  y <- c(1, 0, 1, 0)
  mu <- c(0.8, 0.3, 0.6, 0.2)
  w <- c(1.0, 2.0, 1.5, 0.5)
  expected <- w * (y - mu)
  expect_equal(.compute_residuals("binomial", y, mu, w), expected)
})

test_that("residuals are zero when y equals mu", {
  y <- c(0.3, 0.7)
  mu <- c(0.3, 0.7)
  w <- c(1.0, 2.0)
  expect_equal(.compute_residuals("gaussian", y, mu, w), c(0, 0))
})


# --- working_weights_unweighted ---

test_that("unweighted binomial working weights are mu*(1-mu)", {
  mu <- c(0.2, 0.5, 0.8)
  expected <- mu * (1 - mu)
  expect_equal(.working_weights_unweighted("binomial", mu), expected)
})

test_that("unweighted gaussian working weights are 1/sigma_e^2", {
  mu <- c(1.0, 2.0, 3.0)  # values don't matter for gaussian

  sigma_e <- 2.0
  expected <- rep(1 / 4, 3)
  expect_equal(.working_weights_unweighted("gaussian", mu, sigma_e = sigma_e),
               expected)
})
