# helper-fixtures.R
# Test fixture generators for svyder tests
# ---------------------------------------------------------------------------
# Each fixture returns a self-contained list with all inputs needed for
# DER computation, plus expected analytical values where available.
# No external data files or optional packages required.


#' Balanced Gaussian fixture
#'
#' J=5 groups, n_j=20 each, N=100, p=1 (intercept only).
#' Equal weights, sigma_e=1.0, sigma_theta=0.5.
#' This fixture has known analytical DER values.
#'
#' Analytical derivation for the intercept (beta):
#'   - Working weights: v_i = w_i / sigma_e^2 = 1.0 (equal weights, sigma_e=1)
#'   - H_obs beta-beta block: sum(v_i * x_i^2) + 1/prior_sd^2 = 100 + 1/25 = 100.04
#'   - H_obs beta-theta_j block: sum_{i in j}(v_i * x_i) = 20
#'   - H_obs theta_j-theta_j: sum_{i in j}(v_i) + 1/sigma_theta^2 = 20 + 4 = 24
#'
#' Shrinkage factor B_j:
#'   sigma_theta^2 / (sigma_theta^2 + sigma_e^2/n_j) = 0.25 / (0.25 + 0.05) = 5/6
#'
#' @param seed Random seed for reproducibility.
#' @return Named list with data, draws, and expected analytical values.
make_balanced_gaussian <- function(seed = 42) {

  set.seed(seed)

  # --- Design parameters ---
  J     <- 5L
  n_j   <- 20L
  N     <- J * n_j         # 100
  p     <- 1L              # intercept only
  d     <- p + J           # 6

  sigma_e     <- 1.0
  sigma_theta <- 0.5
  beta_true   <- 2.0
  beta_prior_sd <- 5.0

  # --- Group structure ---
  group <- rep(seq_len(J), each = n_j)

  # --- Random effects ---
  theta_true <- rnorm(J, mean = 0, sd = sigma_theta)

  # --- Design matrix (intercept only) ---
  X <- matrix(1, nrow = N, ncol = p)
  colnames(X) <- "intercept"

  # --- Equal weights ---
  w <- rep(1.0, N)

  # --- PSU = group for this simple case ---
  psu <- group

  # --- Generate continuous outcome ---
  eta <- X %*% beta_true + theta_true[group]
  y   <- as.numeric(eta + rnorm(N, 0, sigma_e))

  # --- Point estimates (use true values for fixture stability) ---
  beta_hat        <- beta_true
  theta_hat       <- theta_true
  sigma_theta_hat <- sigma_theta

  # --- Working weights and residuals for Gaussian ---
  # For Gaussian: working weight = w / sigma_e^2
  # Residual contribution: w * (y - mu) / sigma_e^2 ... but the standalone

  # code uses r = w * (y - mu) with the weight v = w * wt where wt is
  # model-specific. For Gaussian with known sigma_e:
  #   v_i = w_i / sigma_e^2 = 1.0  (working weights for information)
  #   r_i = w_i * (y_i - mu_i) / sigma_e^2  (weighted residuals)
  # But to keep the fixture general, we store raw components.

  mu <- as.numeric(X %*% beta_hat + theta_hat[group])

  # --- Build H_obs analytically ---
  # For Gaussian with sigma_e=1 and w=1:
  #   v = w / sigma_e^2 = 1
  # H_obs[1,1] = sum(v * x^2) + 1/prior_sd^2 = N + 1/25
  # H_obs[1, p+j] = sum_{i in j}(v * x) = n_j
  # H_obs[p+j, p+j] = sum_{i in j}(v) + 1/sigma_theta^2 = n_j + 1/sigma_theta^2

  beta_prior_prec  <- 1 / beta_prior_sd^2
  theta_prior_prec <- 1 / sigma_theta^2

  H_obs <- matrix(0, d, d)
  # beta-beta block
  H_obs[1, 1] <- N + beta_prior_prec
  # beta-theta and theta-beta blocks
  for (j in seq_len(J)) {
    H_obs[1, p + j]     <- n_j
    H_obs[p + j, 1]     <- n_j
    H_obs[p + j, p + j] <- n_j + theta_prior_prec
  }

  # --- Posterior covariance = H_obs^{-1} ---
  H_obs_inv <- solve(H_obs)

  # --- Generate posterior draws via Cholesky decomposition ---
  M <- 4000L
  point_est <- c(beta_hat, theta_hat)

  # Cholesky of posterior covariance
  L <- chol(H_obs_inv)  # upper triangular: H_obs_inv = t(L) %*% L
  Z <- matrix(rnorm(M * d), nrow = M, ncol = d)
  draws_all <- sweep(Z %*% L, 2, point_est, "+")

  draws_beta  <- draws_all[, 1:p, drop = FALSE]
  draws_theta <- draws_all[, (p + 1):d, drop = FALSE]

  colnames(draws_beta)  <- paste0("beta[", seq_len(p), "]")
  colnames(draws_theta) <- paste0("theta[", seq_len(J), "]")

  # --- Parameter types ---
  param_types <- "fe_between"

  # --- Analytical expected values ---

  # Shrinkage factor B_j (identical for all j since balanced):
  #   B = sigma_theta^2 / (sigma_theta^2 + sigma_e^2 / n_j)
  #   B = 0.25 / (0.25 + 1/20) = 0.25 / 0.30 = 5/6
  B_analytical <- sigma_theta^2 / (sigma_theta^2 + sigma_e^2 / n_j)

  # DEFF_j = 1.0 for equal weights (balanced design)
  deff_analytical <- 1.0

  # For the random effects (theta_j), the analytical DER depends on the
  # full sandwich computation. With equal weights and balanced design,
  # V_sand should equal H_obs_inv (since J_cluster = H_obs for the true
  # model), so DER_theta should be approximately 1.0.
  #
  # The exact DER for theta depends on the realized data through J_cluster,
  # so we don't set a hard expected value here -- but we know:
  #   - DER should be close to 1.0 for equal weights (no design effect)
  #   - DER for random effects should be < DER for fixed effects
  #     when there is shrinkage

  list(
    y               = y,
    X               = X,
    group           = group,
    w               = w,
    psu             = psu,
    draws_beta      = draws_beta,
    draws_theta     = draws_theta,
    beta_hat        = beta_hat,
    theta_hat       = theta_hat,
    sigma_theta_hat = sigma_theta_hat,
    sigma_e         = sigma_e,
    param_types     = param_types,
    beta_prior_sd   = beta_prior_sd,
    # Derived quantities
    H_obs           = H_obs,
    H_obs_inv       = H_obs_inv,
    mu              = mu,
    point_est       = point_est,
    # Analytical expectations
    B_analytical    = B_analytical,
    deff_analytical = deff_analytical,
    # Dimensions
    J = J, N = N, p = p, d = d, n_j = n_j, M = M,
    family = "gaussian"
  )
}


#' Unbalanced binomial fixture
#'
#' J=10 groups, N=500, p=2 (intercept + one covariate), binary y.
#' Unequal weights (log-normal distributed).
#' Two covariate types: intercept (between) + continuous (within).
#'
#' @param seed Random seed for reproducibility.
#' @return Named list with data and draws.
make_unbalanced_binomial <- function(seed = 42) {

  set.seed(seed)

  # --- Design parameters ---
  J <- 10L
  N <- 500L
  p <- 2L              # intercept + 1 covariate
  d <- p + J           # 12

  sigma_theta   <- 0.8
  beta_true     <- c(0.5, -0.3)   # intercept, x1
  beta_prior_sd <- 5.0

  # --- Unbalanced group structure ---
  # Varying cluster sizes: 30-70 per group
  n_per_group <- rmultinom(1, N, prob = rep(1, J))[, 1]
  # Ensure no empty groups
  n_per_group[n_per_group == 0] <- 1L
  # Adjust to get exactly N
  diff_n <- N - sum(n_per_group)
  if (diff_n > 0) {
    add_to <- sample(J, diff_n, replace = TRUE)
    for (idx in add_to) n_per_group[idx] <- n_per_group[idx] + 1L
  } else if (diff_n < 0) {
    sub_from <- sample(which(n_per_group > 1), abs(diff_n), replace = TRUE)
    for (idx in sub_from) n_per_group[idx] <- n_per_group[idx] - 1L
  }

  group <- rep(seq_len(J), times = n_per_group)

  # --- Random effects ---
  theta_true <- rnorm(J, mean = 0, sd = sigma_theta)

  # --- Design matrix ---
  x1 <- rnorm(N, mean = 0, sd = 1)
  X  <- cbind(intercept = 1, x1 = x1)

  # --- Unequal weights (log-normal) ---
  w <- exp(rnorm(N, mean = 0, sd = 0.5))
  w <- w / mean(w)   # normalize to mean 1

  # --- PSU = group ---
  psu <- group

  # --- Binary outcome ---
  eta   <- X %*% beta_true + theta_true[group]
  prob  <- 1 / (1 + exp(-eta))
  y     <- rbinom(N, size = 1, prob = prob)

  # --- Point estimates (use true values) ---
  beta_hat        <- beta_true
  theta_hat       <- theta_true
  sigma_theta_hat <- sigma_theta

  # --- Build H_obs for logistic model ---
  mu <- as.numeric(1 / (1 + exp(-(X %*% beta_hat + theta_hat[group]))))
  wt <- mu * (1 - mu)
  v  <- w * wt

  beta_prior_prec  <- 1 / beta_prior_sd^2
  theta_prior_prec <- 1 / sigma_theta^2

  H_obs <- matrix(0, d, d)
  H_obs[1:p, 1:p] <- crossprod(X, X * v)
  for (j in seq_len(J)) {
    idx_j <- which(group == j)
    v_j   <- v[idx_j]
    if (length(idx_j) == 1L) {
      bt_j <- v_j * X[idx_j, ]
    } else {
      bt_j <- colSums(X[idx_j, , drop = FALSE] * v_j)
    }
    H_obs[1:p, p + j]     <- bt_j
    H_obs[p + j, 1:p]     <- bt_j
    H_obs[p + j, p + j]   <- sum(v_j)
  }
  for (k in seq_len(p)) {
    H_obs[k, k] <- H_obs[k, k] + beta_prior_prec
  }
  for (j in seq_len(J)) {
    H_obs[p + j, p + j] <- H_obs[p + j, p + j] + theta_prior_prec
  }

  # --- Posterior covariance approximation (Laplace) ---
  H_obs_inv <- solve(H_obs)

  # --- Generate posterior draws via Cholesky ---
  M <- 4000L
  point_est <- c(beta_hat, theta_hat)
  L <- chol(H_obs_inv)
  Z <- matrix(rnorm(M * d), nrow = M, ncol = d)
  draws_all <- sweep(Z %*% L, 2, point_est, "+")

  draws_beta  <- draws_all[, 1:p, drop = FALSE]
  draws_theta <- draws_all[, (p + 1):d, drop = FALSE]

  colnames(draws_beta)  <- paste0("beta[", seq_len(p), "]")
  colnames(draws_theta) <- paste0("theta[", seq_len(J), "]")

  # --- Parameter types ---
  param_types <- c("fe_between", "fe_within")

  list(
    y               = y,
    X               = X,
    group           = group,
    w               = w,
    psu             = psu,
    draws_beta      = draws_beta,
    draws_theta     = draws_theta,
    beta_hat        = beta_hat,
    theta_hat       = theta_hat,
    sigma_theta_hat = sigma_theta_hat,
    param_types     = param_types,
    beta_prior_sd   = beta_prior_sd,
    H_obs           = H_obs,
    H_obs_inv       = H_obs_inv,
    mu              = mu,
    v               = v,
    wt              = wt,
    point_est       = point_est,
    J = J, N = N, p = p, d = d, M = M,
    family = "binomial"
  )
}


#' Minimal J=2 fixture
#'
#' J=2 groups, n_j=10, N=20, p=1 (intercept only).
#' Minimum viable cluster count to test edge cases.
#'
#' @param seed Random seed for reproducibility.
#' @return Named list with data and draws.
make_minimal_j2 <- function(seed = 42) {

  set.seed(seed)

  # --- Design parameters ---
  J   <- 2L
  n_j <- 10L
  N   <- J * n_j        # 20
  p   <- 1L             # intercept only
  d   <- p + J          # 3

  sigma_theta   <- 0.3
  beta_true     <- 1.0
  sigma_e       <- 1.0
  beta_prior_sd <- 5.0

  # --- Group structure ---
  group <- rep(seq_len(J), each = n_j)

  # --- Random effects ---
  theta_true <- rnorm(J, mean = 0, sd = sigma_theta)

  # --- Design matrix ---
  X <- matrix(1, nrow = N, ncol = p)
  colnames(X) <- "intercept"

  # --- Equal weights ---
  w <- rep(1.0, N)

  # --- PSU = group ---
  psu <- group

  # --- Generate continuous outcome ---
  eta <- X %*% beta_true + theta_true[group]
  y   <- as.numeric(eta + rnorm(N, 0, sigma_e))

  # --- Point estimates ---
  beta_hat        <- beta_true
  theta_hat       <- theta_true
  sigma_theta_hat <- sigma_theta

  mu <- as.numeric(X %*% beta_hat + theta_hat[group])

  # --- Build H_obs (Gaussian, sigma_e = 1, w = 1) ---
  beta_prior_prec  <- 1 / beta_prior_sd^2
  theta_prior_prec <- 1 / sigma_theta^2

  H_obs <- matrix(0, d, d)
  H_obs[1, 1] <- N + beta_prior_prec
  for (j in seq_len(J)) {
    H_obs[1, p + j]     <- n_j
    H_obs[p + j, 1]     <- n_j
    H_obs[p + j, p + j] <- n_j + theta_prior_prec
  }

  H_obs_inv <- solve(H_obs)

  # --- Generate posterior draws via Cholesky ---
  M <- 2000L
  point_est <- c(beta_hat, theta_hat)
  L <- chol(H_obs_inv)
  Z <- matrix(rnorm(M * d), nrow = M, ncol = d)
  draws_all <- sweep(Z %*% L, 2, point_est, "+")

  draws_beta  <- draws_all[, 1:p, drop = FALSE]
  draws_theta <- draws_all[, (p + 1):d, drop = FALSE]

  colnames(draws_beta)  <- paste0("beta[", seq_len(p), "]")
  colnames(draws_theta) <- paste0("theta[", seq_len(J), "]")

  # --- Parameter types ---
  param_types <- "fe_between"

  list(
    y               = y,
    X               = X,
    group           = group,
    w               = w,
    psu             = psu,
    draws_beta      = draws_beta,
    draws_theta     = draws_theta,
    beta_hat        = beta_hat,
    theta_hat       = theta_hat,
    sigma_theta_hat = sigma_theta_hat,
    sigma_e         = sigma_e,
    param_types     = param_types,
    beta_prior_sd   = beta_prior_sd,
    H_obs           = H_obs,
    H_obs_inv       = H_obs_inv,
    mu              = mu,
    point_est       = point_est,
    J = J, N = N, p = p, d = d, n_j = n_j, M = M,
    family = "gaussian"
  )
}
