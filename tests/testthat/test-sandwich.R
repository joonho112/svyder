# test-sandwich.R
# Tests for sandwich variance estimation internals
# (.build_H_obs, .safe_invert, .build_J_cluster, .build_V_sand)


# ============================================================================
# .build_H_obs tests
# ============================================================================

test_that("H_obs is symmetric for balanced gaussian", {
  fix <- make_balanced_gaussian()

  # Working weights for gaussian with sigma_e=1, w=1: v = 1
  v <- rep(1.0, fix$N)

  H <- svyder:::.build_H_obs(
    X = fix$X, v = v, group = fix$group,
    J = fix$J, p = fix$p,
    beta_prior_prec  = 1 / fix$beta_prior_sd^2,
    theta_prior_prec = 1 / fix$sigma_theta_hat^2
  )

  expect_equal(H, t(H), tolerance = 1e-12,
               info = "H_obs must be symmetric")
})

test_that("H_obs is symmetric for unbalanced binomial", {
  fix <- make_unbalanced_binomial()

  H <- svyder:::.build_H_obs(
    X = fix$X, v = fix$v, group = fix$group,
    J = fix$J, p = fix$p,
    beta_prior_prec  = 1 / fix$beta_prior_sd^2,
    theta_prior_prec = 1 / fix$sigma_theta_hat^2
  )

  expect_equal(H, t(H), tolerance = 1e-12)
})

test_that("H_obs is positive definite for balanced gaussian", {
  fix <- make_balanced_gaussian()
  v <- rep(1.0, fix$N)

  H <- svyder:::.build_H_obs(
    X = fix$X, v = v, group = fix$group,
    J = fix$J, p = fix$p,
    beta_prior_prec  = 1 / fix$beta_prior_sd^2,
    theta_prior_prec = 1 / fix$sigma_theta_hat^2
  )

  eigenvalues <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigenvalues > 0),
              info = "All eigenvalues of H_obs must be positive")
})

test_that("H_obs is positive definite for unbalanced binomial", {
  fix <- make_unbalanced_binomial()

  H <- svyder:::.build_H_obs(
    X = fix$X, v = fix$v, group = fix$group,
    J = fix$J, p = fix$p,
    beta_prior_prec  = 1 / fix$beta_prior_sd^2,
    theta_prior_prec = 1 / fix$sigma_theta_hat^2
  )

  eigenvalues <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigenvalues > 0))
})

test_that("H_obs matches analytical values for balanced gaussian", {
  fix <- make_balanced_gaussian()
  v <- rep(1.0, fix$N)

  H <- svyder:::.build_H_obs(
    X = fix$X, v = v, group = fix$group,
    J = fix$J, p = fix$p,
    beta_prior_prec  = 1 / fix$beta_prior_sd^2,
    theta_prior_prec = 1 / fix$sigma_theta_hat^2
  )

  # Check against the fixture's pre-computed H_obs
  expect_equal(H, fix$H_obs, tolerance = 1e-10,
               info = "H_obs from .build_H_obs must match analytical H_obs")
})

test_that("H_obs matches standalone code for unbalanced binomial", {
  fix <- make_unbalanced_binomial()

  H <- svyder:::.build_H_obs(
    X = fix$X, v = fix$v, group = fix$group,
    J = fix$J, p = fix$p,
    beta_prior_prec  = 1 / fix$beta_prior_sd^2,
    theta_prior_prec = 1 / fix$sigma_theta_hat^2
  )

  expect_equal(H, fix$H_obs, tolerance = 1e-10)
})

test_that("H_obs handles zero prior precision (flat prior)", {
  fix <- make_balanced_gaussian()
  v <- rep(1.0, fix$N)

  H_flat <- svyder:::.build_H_obs(
    X = fix$X, v = v, group = fix$group,
    J = fix$J, p = fix$p,
    beta_prior_prec  = 0,
    theta_prior_prec = 1 / fix$sigma_theta_hat^2
  )

  # Beta-beta block should NOT include prior precision
  expect_equal(H_flat[1, 1], fix$N,
               info = "Flat prior means H_obs[1,1] = N (no prior term)")
})

test_that("known-answer test for small case", {
  # p=1, J=3, n_j=5, intercept only, equal weights, v=1
  J <- 3L
  p <- 1L
  n_j <- 5L
  N <- J * n_j  # 15
  d <- p + J    # 4

  X <- matrix(1, nrow = N, ncol = p)
  v <- rep(1.0, N)
  group <- rep(seq_len(J), each = n_j)

  beta_prior_prec  <- 0.04   # 1/5^2
  theta_prior_prec <- 4.0    # 1/0.5^2

  H <- svyder:::.build_H_obs(X, v, group, J, p, beta_prior_prec, theta_prior_prec)

  # Hand computation:
  # H[1,1] = sum(v * x^2) + beta_prior_prec = 15 + 0.04 = 15.04
  expect_equal(H[1, 1], 15.04, tolerance = 1e-12)

  # H[1, p+j] = sum_{i in j}(v * x) = 5 for all j
  for (j in 1:3) {
    expect_equal(H[1, p + j], 5.0, tolerance = 1e-12)
    expect_equal(H[p + j, 1], 5.0, tolerance = 1e-12)
  }

  # H[p+j, p+j] = sum_{i in j}(v) + theta_prior_prec = 5 + 4 = 9
  for (j in 1:3) {
    expect_equal(H[p + j, p + j], 9.0, tolerance = 1e-12)
  }

  # Off-diagonal theta blocks should be 0
  expect_equal(H[p + 1, p + 2], 0.0, tolerance = 1e-12)
  expect_equal(H[p + 2, p + 3], 0.0, tolerance = 1e-12)
  expect_equal(H[p + 1, p + 3], 0.0, tolerance = 1e-12)
})


# ============================================================================
# .build_J_cluster tests
# ============================================================================

test_that("J_cluster is symmetric for balanced gaussian", {
  fix <- make_balanced_gaussian()

  # Weighted residuals for gaussian: r = w * (y - mu) / sigma_e^2
  r <- fix$w * (fix$y - fix$mu) / fix$sigma_e^2

  Jc <- svyder:::.build_J_cluster(
    X = fix$X, r = r, psu = fix$psu, group = fix$group,
    p = fix$p, J = fix$J
  )

  expect_equal(Jc, t(Jc), tolerance = 1e-12,
               info = "J_cluster must be symmetric")
})

test_that("J_cluster is symmetric for unbalanced binomial", {
  fix <- make_unbalanced_binomial()

  # Weighted residuals for logistic: r = w * (y - mu)
  r <- fix$w * (fix$y - fix$mu)

  Jc <- svyder:::.build_J_cluster(
    X = fix$X, r = r, psu = fix$psu, group = fix$group,
    p = fix$p, J = fix$J
  )

  expect_equal(Jc, t(Jc), tolerance = 1e-12)
})

test_that("J_cluster is positive semi-definite", {
  fix <- make_balanced_gaussian()
  r <- fix$w * (fix$y - fix$mu) / fix$sigma_e^2

  Jc <- svyder:::.build_J_cluster(
    X = fix$X, r = r, psu = fix$psu, group = fix$group,
    p = fix$p, J = fix$J
  )

  eigenvalues <- eigen(Jc, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigenvalues >= -1e-10),
              info = "J_cluster eigenvalues must be non-negative (PSD)")
})

test_that("J_cluster handles single-observation PSUs", {
  # Each observation is its own PSU
  set.seed(99)
  J <- 2L; p <- 1L; N <- 10L
  X <- matrix(1, N, p)
  r <- rnorm(N)
  group <- rep(seq_len(J), each = N / J)
  psu <- seq_len(N)  # each obs is its own PSU

  Jc <- svyder:::.build_J_cluster(X, r, psu, group, p, J)

  expect_equal(Jc, t(Jc), tolerance = 1e-12)
  # With single-obs PSUs, J_cluster should equal sum of r_i^2 x_i x_i^T
  # for the beta block
  expect_equal(Jc[1, 1], sum(r^2), tolerance = 1e-10,
               info = "Single-obs PSUs: J_cluster[1,1] = sum(r^2)")
})


# ============================================================================
# .build_V_sand tests
# ============================================================================

test_that("V_sand = H_inv J H_inv identity", {
  fix <- make_balanced_gaussian()

  v <- rep(1.0, fix$N)
  r <- fix$w * (fix$y - fix$mu) / fix$sigma_e^2

  H <- svyder:::.build_H_obs(
    X = fix$X, v = v, group = fix$group,
    J = fix$J, p = fix$p,
    beta_prior_prec  = 1 / fix$beta_prior_sd^2,
    theta_prior_prec = 1 / fix$sigma_theta_hat^2
  )

  H_inv <- solve(H)

  Jc <- svyder:::.build_J_cluster(
    X = fix$X, r = r, psu = fix$psu, group = fix$group,
    p = fix$p, J = fix$J
  )

  V <- svyder:::.build_V_sand(H_inv, Jc)

  # Verify the sandwich formula
  V_expected <- H_inv %*% Jc %*% H_inv
  V_expected <- (V_expected + t(V_expected)) / 2

  expect_equal(V, V_expected, tolerance = 1e-12,
               info = "V_sand must equal H_inv %*% J %*% H_inv (symmetrized)")
})

test_that("V_sand is symmetric", {
  fix <- make_unbalanced_binomial()

  r <- fix$w * (fix$y - fix$mu)

  H <- svyder:::.build_H_obs(
    X = fix$X, v = fix$v, group = fix$group,
    J = fix$J, p = fix$p,
    beta_prior_prec  = 1 / fix$beta_prior_sd^2,
    theta_prior_prec = 1 / fix$sigma_theta_hat^2
  )

  H_inv <- solve(H)

  Jc <- svyder:::.build_J_cluster(
    X = fix$X, r = r, psu = fix$psu, group = fix$group,
    p = fix$p, J = fix$J
  )

  V <- svyder:::.build_V_sand(H_inv, Jc)

  expect_equal(V, t(V), tolerance = 1e-12)
})

test_that("V_sand has non-negative diagonal", {
  fix <- make_balanced_gaussian()
  v <- rep(1.0, fix$N)
  r <- fix$w * (fix$y - fix$mu) / fix$sigma_e^2

  H <- svyder:::.build_H_obs(
    X = fix$X, v = v, group = fix$group,
    J = fix$J, p = fix$p,
    beta_prior_prec  = 1 / fix$beta_prior_sd^2,
    theta_prior_prec = 1 / fix$sigma_theta_hat^2
  )
  H_inv <- solve(H)

  Jc <- svyder:::.build_J_cluster(
    X = fix$X, r = r, psu = fix$psu, group = fix$group,
    p = fix$p, J = fix$J
  )

  V <- svyder:::.build_V_sand(H_inv, Jc)

  expect_true(all(diag(V) >= -1e-15),
              info = "Diagonal elements of V_sand must be non-negative")
})


# ============================================================================
# .safe_invert tests
# ============================================================================

test_that("safe_invert inverts a well-conditioned matrix", {
  A <- matrix(c(2, 1, 1, 3), 2, 2)
  A_inv <- svyder:::.safe_invert(A)

  expect_equal(A %*% A_inv, diag(2), tolerance = 1e-12,
               info = "A %*% A_inv should be identity")
})

test_that("safe_invert handles near-singular matrix", {
  skip_if_not_installed("Matrix")

  # Create a near-singular matrix
  d <- 4
  A <- diag(d)
  A[1, 1] <- 1e-16  # nearly singular

  # Should not error; falls back to nearPD
  A_inv <- svyder:::.safe_invert(A)

  expect_true(is.matrix(A_inv))
  expect_equal(nrow(A_inv), d)
  expect_equal(ncol(A_inv), d)
})

test_that("safe_invert returns correct result for identity matrix", {
  d <- 5
  I_inv <- svyder:::.safe_invert(diag(d))
  expect_equal(I_inv, diag(d), tolerance = 1e-12)
})

test_that("safe_invert returns correct result for diagonal matrix", {
  d <- 4
  diag_vals <- c(2, 5, 0.5, 10)
  A <- diag(diag_vals)
  A_inv <- svyder:::.safe_invert(A)

  expect_equal(diag(A_inv), 1 / diag_vals, tolerance = 1e-12)
})


# ============================================================================
# Full pipeline integration test
# ============================================================================

test_that("full sandwich pipeline matches standalone code output", {
  # Use the unbalanced binomial fixture to test exact match
  # with the standalone compute_der.R logic
  fix <- make_unbalanced_binomial()

  v <- fix$v  # pre-computed: w * mu * (1 - mu)
  r <- fix$w * (fix$y - fix$mu)

  beta_prior_prec  <- 1 / fix$beta_prior_sd^2
  theta_prior_prec <- 1 / fix$sigma_theta_hat^2

  # Build via modular functions
  H <- svyder:::.build_H_obs(
    X = fix$X, v = v, group = fix$group,
    J = fix$J, p = fix$p,
    beta_prior_prec  = beta_prior_prec,
    theta_prior_prec = theta_prior_prec
  )

  H_inv <- svyder:::.safe_invert(H)

  Jc <- svyder:::.build_J_cluster(
    X = fix$X, r = r, psu = fix$psu, group = fix$group,
    p = fix$p, J = fix$J
  )

  V <- svyder:::.build_V_sand(H_inv, Jc)

  # Replicate standalone code inline for comparison
  d <- fix$d
  p <- fix$p
  J <- fix$J
  N <- fix$N

  H_standalone <- matrix(0, d, d)
  H_standalone[1:p, 1:p] <- crossprod(fix$X, fix$X * v)
  for (j in seq_len(J)) {
    idx_j <- which(fix$group == j)
    v_j <- v[idx_j]
    if (length(idx_j) == 1L) {
      bt_j <- v_j * fix$X[idx_j, ]
    } else {
      bt_j <- colSums(fix$X[idx_j, , drop = FALSE] * v_j)
    }
    H_standalone[1:p, p + j] <- bt_j
    H_standalone[p + j, 1:p] <- bt_j
    H_standalone[p + j, p + j] <- sum(v_j)
  }
  for (k in seq_len(p)) {
    H_standalone[k, k] <- H_standalone[k, k] + beta_prior_prec
  }
  for (j in seq_len(J)) {
    H_standalone[p + j, p + j] <- H_standalone[p + j, p + j] + theta_prior_prec
  }

  H_standalone_inv <- solve(H_standalone)

  G <- max(fix$psu)
  J_standalone <- matrix(0, d, d)
  for (g in seq_len(G)) {
    idx_g <- which(fix$psu == g)
    if (length(idx_g) == 0L) next
    s_g <- numeric(d)
    if (length(idx_g) == 1L) {
      s_g[1:p] <- r[idx_g] * fix$X[idx_g, ]
    } else {
      s_g[1:p] <- colSums(fix$X[idx_g, , drop = FALSE] * r[idx_g])
    }
    for (j_state in unique(fix$group[idx_g])) {
      idx_gj <- idx_g[fix$group[idx_g] == j_state]
      s_g[p + j_state] <- s_g[p + j_state] + sum(r[idx_gj])
    }
    J_standalone <- J_standalone + tcrossprod(s_g)
  }
  J_standalone <- (J_standalone + t(J_standalone)) / 2

  V_standalone <- H_standalone_inv %*% J_standalone %*% H_standalone_inv
  V_standalone <- (V_standalone + t(V_standalone)) / 2

  # Compare
  expect_equal(H, H_standalone, tolerance = 1e-12,
               info = "H_obs must match standalone code exactly")
  expect_equal(H_inv, H_standalone_inv, tolerance = 1e-10,
               info = "H_obs_inv must match standalone code")
  expect_equal(Jc, J_standalone, tolerance = 1e-12,
               info = "J_cluster must match standalone code exactly")
  expect_equal(V, V_standalone, tolerance = 1e-10,
               info = "V_sand must match standalone code")
})

test_that("full pipeline works for minimal J=2 fixture", {
  fix <- make_minimal_j2()

  v <- rep(1.0, fix$N)  # gaussian with sigma_e=1, w=1
  r <- fix$w * (fix$y - fix$mu) / fix$sigma_e^2

  H <- svyder:::.build_H_obs(
    X = fix$X, v = v, group = fix$group,
    J = fix$J, p = fix$p,
    beta_prior_prec  = 1 / fix$beta_prior_sd^2,
    theta_prior_prec = 1 / fix$sigma_theta_hat^2
  )

  expect_equal(H, fix$H_obs, tolerance = 1e-10)

  H_inv <- svyder:::.safe_invert(H)
  expect_equal(H %*% H_inv, diag(fix$d), tolerance = 1e-10)

  Jc <- svyder:::.build_J_cluster(
    X = fix$X, r = r, psu = fix$psu, group = fix$group,
    p = fix$p, J = fix$J
  )

  expect_equal(Jc, t(Jc), tolerance = 1e-12)

  V <- svyder:::.build_V_sand(H_inv, Jc)

  expect_equal(V, t(V), tolerance = 1e-12)
  expect_true(all(diag(V) >= -1e-15))
})
