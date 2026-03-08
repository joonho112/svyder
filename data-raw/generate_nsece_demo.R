###############################################################################
# generate_nsece_demo.R
# Generate the nsece_demo bundled dataset for the svyder package
# ---------------------------------------------------------------------------
# Ports the standalone generate_nsece_demo() function and creates synthetic
# posterior draws from a multivariate normal approximation.
###############################################################################

# --- Generate synthetic NSECE-like survey data ---

generate_nsece_demo_data <- function(seed = 42) {
  set.seed(seed)

  # Actual NSECE 2019 state sample sizes (J=51 states including DC)
  state_n <- c(17, 18, 22, 24, 26, 27, 27, 34, 36, 38, 40,
               40, 43, 43, 46, 50, 50, 51, 52, 56, 57, 57, 58, 59, 61,
               64, 74, 75, 77, 86, 98, 102, 103, 110, 117, 120, 121,
               131, 138, 141, 177, 178, 183, 186, 192, 201, 362, 471,
               558, 578, 1110)

  N <- sum(state_n)
  J <- length(state_n)

  # Group indicator
  group <- rep(seq_len(J), times = state_n)

  # Covariates
  poverty_raw <- rnorm(N, mean = 0, sd = 1)
  poverty_cwc <- numeric(N)
  for (j in seq_len(J)) {
    idx_j <- which(group == j)
    poverty_cwc[idx_j] <- poverty_raw[idx_j] - mean(poverty_raw[idx_j])
  }

  tiered_reim <- rep(0, J)
  tiered_reim[sample(J, size = round(J * 0.55))] <- 1
  x_between <- tiered_reim[group]

  X <- cbind(1, poverty_cwc, x_between)
  colnames(X) <- c("intercept", "poverty_cwc", "tiered_reim")

  # True parameters
  beta_true <- c(0.25, -0.15, 0.16)
  sigma_theta <- 0.66

  # Random effects
  theta_true <- rnorm(J, mean = 0, sd = sigma_theta)

  # Linear predictor and outcome
  eta <- as.numeric(X %*% beta_true) + theta_true[group]
  mu <- 1 / (1 + exp(-eta))
  y <- rbinom(N, size = 1, prob = mu)

  # Survey weights (log-normal, normalized within state)
  sigma_w <- sqrt(log(3))
  log_w <- rnorm(N, mean = 0, sd = sigma_w)
  w_raw <- exp(log_w)
  w <- numeric(N)
  for (j in seq_len(J)) {
    idx_j <- which(group == j)
    n_j <- length(idx_j)
    w[idx_j] <- w_raw[idx_j] / sum(w_raw[idx_j]) * n_j
  }

  # PSU structure (8-16 observations per PSU)
  psu <- integer(N)
  psu_counter <- 0L
  for (j in seq_len(J)) {
    idx_j <- which(group == j)
    n_j <- length(idx_j)
    psu_size <- sample(8:16, size = 1)
    n_psus_j <- max(1L, ceiling(n_j / psu_size))
    psu_labels <- rep(seq_len(n_psus_j), each = psu_size, length.out = n_j)
    psu[idx_j] <- psu_counter + psu_labels
    psu_counter <- psu_counter + n_psus_j
  }

  list(
    y = y, X = X, group = group, w = w, psu = psu,
    beta_true = beta_true, sigma_theta = sigma_theta,
    theta_true = theta_true, state_n = state_n,
    N = N, J = J, p = 3L
  )
}


# --- Build posterior draws from Laplace approximation ---

build_posterior_draws <- function(data, S = 4000L, beta_prior_sd = 5.0) {

  y     <- data$y
  X     <- data$X
  group <- data$group
  w     <- data$w
  N     <- data$N
  J     <- data$J
  p     <- data$p
  d     <- p + J

  beta_true   <- data$beta_true
  theta_true  <- data$theta_true
  sigma_theta <- data$sigma_theta

  # Point estimates (use true values for stability)
  beta_hat  <- beta_true
  theta_hat <- theta_true

  # Fitted values
  eta <- as.numeric(X %*% beta_hat) + theta_hat[group]
  mu  <- 1 / (1 + exp(-eta))

  # Working weights (binomial): v_i = w_i * mu_i * (1 - mu_i)
  wt <- mu * (1 - mu)
  v  <- w * wt

  # Prior precisions
  beta_prior_prec  <- 1 / beta_prior_sd^2
  theta_prior_prec <- 1 / sigma_theta^2

  # Build observed information matrix
  H_obs <- matrix(0, d, d)
  H_obs[1:p, 1:p] <- crossprod(X, X * v)

  for (j in seq_len(J)) {
    idx_j <- which(group == j)
    v_j <- v[idx_j]
    if (length(idx_j) == 1L) {
      bt_j <- v_j * X[idx_j, ]
    } else {
      bt_j <- colSums(X[idx_j, , drop = FALSE] * v_j)
    }
    H_obs[1:p, p + j] <- bt_j
    H_obs[p + j, 1:p] <- bt_j
    H_obs[p + j, p + j] <- sum(v_j)
  }

  for (k in seq_len(p)) {
    H_obs[k, k] <- H_obs[k, k] + beta_prior_prec
  }
  for (j in seq_len(J)) {
    H_obs[p + j, p + j] <- H_obs[p + j, p + j] + theta_prior_prec
  }

  # Posterior covariance (Laplace approximation)
  H_obs_inv <- solve(H_obs)

  # Generate draws via Cholesky decomposition
  point_est <- c(beta_hat, theta_hat)
  L <- chol(H_obs_inv)  # upper Cholesky: H_obs_inv = t(L) %*% L
  Z <- matrix(rnorm(S * d), nrow = S, ncol = d)
  draws <- sweep(Z %*% L, 2, point_est, "+")

  colnames(draws) <- c(paste0("beta[", seq_len(p), "]"),
                        paste0("theta[", seq_len(J), "]"))

  draws
}


# --- Main: generate and save ---

cat("Generating nsece_demo dataset...\n")

raw_data <- generate_nsece_demo_data(seed = 42)
draws    <- build_posterior_draws(raw_data, S = 4000L)

nsece_demo <- list(
  draws       = draws,
  y           = raw_data$y,
  X           = raw_data$X,
  group       = raw_data$group,
  weights     = raw_data$w,
  psu         = raw_data$psu,
  param_types = c("fe_between", "fe_within", "fe_between"),
  family      = "binomial",
  sigma_theta = raw_data$sigma_theta,
  N           = raw_data$N,
  J           = raw_data$J,
  p           = raw_data$p
)

# Ensure data directory exists
if (!dir.exists("data")) dir.create("data")

save(nsece_demo, file = "data/nsece_demo.rda", compress = "xz")
cat(sprintf("  Saved data/nsece_demo.rda\n"))
cat(sprintf("  N = %d, J = %d, p = %d\n", nsece_demo$N, nsece_demo$J, nsece_demo$p))
cat(sprintf("  draws: %d x %d\n", nrow(nsece_demo$draws), ncol(nsece_demo$draws)))
cat("Done.\n")
