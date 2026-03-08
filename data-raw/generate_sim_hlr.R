###############################################################################
# generate_sim_hlr.R
# Generate the sim_hlr bundled dataset for the svyder package
# ---------------------------------------------------------------------------
# A small Gaussian hierarchical linear regression dataset for quick testing.
# J=10, n_j=20, N=200, p=2 (intercept + 1 within covariate).
# Equal weights => DEFF=1 => DER should be close to 1.0.
###############################################################################

cat("Generating sim_hlr dataset...\n")

set.seed(123)

# --- Design parameters ---
J     <- 10L
n_j   <- 20L
N     <- J * n_j         # 200
p     <- 2L              # intercept + 1 within covariate
d     <- p + J           # 12

sigma_e       <- 1.0
sigma_theta   <- 0.5
beta_true     <- c(2.0, -0.5)   # intercept, x_within
beta_prior_sd <- 5.0

# --- Group structure ---
group <- rep(seq_len(J), each = n_j)

# --- Random effects ---
theta_true <- rnorm(J, mean = 0, sd = sigma_theta)

# --- Design matrix ---
# x_within: group-mean centered covariate (pure within-cluster variation)
x_raw <- rnorm(N)
x_within <- numeric(N)
for (j in seq_len(J)) {
  idx_j <- which(group == j)
  x_within[idx_j] <- x_raw[idx_j] - mean(x_raw[idx_j])
}
X <- cbind(intercept = 1, x_within = x_within)

# --- Equal weights (DEFF = 1) ---
w <- rep(1.0, N)

# --- PSU = group (one PSU per cluster) ---
psu <- group

# --- Generate continuous outcome ---
eta <- as.numeric(X %*% beta_true) + theta_true[group]
y   <- eta + rnorm(N, 0, sigma_e)

# --- Build H_obs analytically (Gaussian, sigma_e = 1, w = 1) ---
# Working weights: v_i = w_i / sigma_e^2 = 1
v <- w / sigma_e^2

beta_prior_prec  <- 1 / beta_prior_sd^2
theta_prior_prec <- 1 / sigma_theta^2

H_obs <- matrix(0, d, d)
H_obs[1:p, 1:p] <- crossprod(X, X * v)

for (j in seq_len(J)) {
  idx_j <- which(group == j)
  v_j <- v[idx_j]
  bt_j <- colSums(X[idx_j, , drop = FALSE] * v_j)
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

# --- Posterior covariance (Laplace approximation) ---
H_obs_inv <- solve(H_obs)

# --- Generate posterior draws ---
S <- 4000L
beta_hat  <- beta_true
theta_hat <- theta_true
point_est <- c(beta_hat, theta_hat)

L <- chol(H_obs_inv)
Z <- matrix(rnorm(S * d), nrow = S, ncol = d)
draws <- sweep(Z %*% L, 2, point_est, "+")

colnames(draws) <- c(paste0("beta[", seq_len(p), "]"),
                      paste0("theta[", seq_len(J), "]"))

# --- Pre-computed reference values ---
# With equal weights, DEFF = 1, so DER should be close to 1.0
# Shrinkage factor: B = sigma_theta^2 / (sigma_theta^2 + sigma_e^2 / n_j)
#                     = 0.25 / (0.25 + 0.05) = 5/6 ~ 0.833
B_ref <- sigma_theta^2 / (sigma_theta^2 + sigma_e^2 / n_j)

# --- Package into sim_hlr ---
sim_hlr <- list(
  draws       = draws,
  y           = y,
  X           = X,
  group       = group,
  weights     = w,
  psu         = psu,
  param_types = c("fe_between", "fe_within"),
  family      = "gaussian",
  sigma_theta = sigma_theta,
  sigma_e     = sigma_e,
  N           = N,
  J           = J,
  p           = p,
  # Reference values for testing
  B_ref       = B_ref,
  deff_ref    = 1.0
)

# Ensure data directory exists
if (!dir.exists("data")) dir.create("data")

save(sim_hlr, file = "data/sim_hlr.rda", compress = "xz")
cat(sprintf("  Saved data/sim_hlr.rda\n"))
cat(sprintf("  N = %d, J = %d, p = %d\n", sim_hlr$N, sim_hlr$J, sim_hlr$p))
cat(sprintf("  draws: %d x %d\n", nrow(sim_hlr$draws), ncol(sim_hlr$draws)))
cat(sprintf("  B_ref = %.4f, deff_ref = %.1f\n", sim_hlr$B_ref, sim_hlr$deff_ref))
cat("Done.\n")
