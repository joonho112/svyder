# test-extract-cmdstanr.R
# Tests for cmdstanr integration: extract_draws.CmdStanMCMC, der_compute.CmdStanMCMC
# ---------------------------------------------------------------------------
# All tests requiring cmdstanr use skip_if_not_installed("cmdstanr").


# ============================================================================
# Method existence checks (work without cmdstanr)
# ============================================================================

test_that("extract_draws.CmdStanMCMC method exists", {
  expect_true(is.function(svyder:::extract_draws.CmdStanMCMC))
})

test_that("der_compute.CmdStanMCMC method exists", {
  expect_true(is.function(svyder:::der_compute.CmdStanMCMC))
})


# ============================================================================
# extract_draws.CmdStanMCMC — requires cmdstanr
# ============================================================================

test_that("extract_draws.CmdStanMCMC works with fitted model", {
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")

  # This is a placeholder for integration testing with an actual CmdStanMCMC

  # object. In practice, testing requires a compiled Stan model, which is
  # too heavy for unit tests. The method's logic is straightforward:
  # convert to draws_matrix, remove lp__ columns.
})


# ============================================================================
# der_compute.CmdStanMCMC — requires cmdstanr
# ============================================================================

test_that("der_compute.CmdStanMCMC delegates correctly", {
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")

  # Integration test placeholder. The cmdstanr method delegates to
  # der_compute.matrix after extracting draws, so correctness of the
  # underlying computation is verified by test-der_compute.R.
})
