# sandwich.R
# Sandwich variance estimation internals
# ---------------------------------------------------------------------------
# These functions build the three matrices needed for the sandwich estimator:
#   H_obs     : observed information (d x d)
#   J_cluster : clustered score outer product (d x d)
#   V_sand    : sandwich variance H_inv %*% J %*% H_inv (d x d)
#
# EXACT port of compute_der.R lines 64-120, refactored into modular functions.
# The math is identical; only the interface is changed so that family-specific
# quantities (working weights v, weighted residuals r) are passed in rather
# than computed here.
#
# All functions are internal (prefixed with .) and not exported.


#' Build the observed information matrix
#'
#' Constructs the (p+J) x (p+J) block matrix H_obs from working weights.
#'
#' @param X Design matrix (N x p).
#' @param v Working weights (length N).
#' @param group Integer group indicator (1:J), length N.
#' @param J Number of groups.
#' @param p Number of fixed-effect parameters.
#' @param beta_prior_prec Scalar prior precision for beta.
#' @param theta_prior_prec Scalar prior precision for theta.
#'
#' @return Symmetric (p+J) x (p+J) matrix.
#'
#' @keywords internal
.build_H_obs <- function(X, v, group, J, p, beta_prior_prec, theta_prior_prec) {

  d <- p + J
  H_obs <- matrix(0, d, d)

  # --- Beta-beta block: X^T diag(v) X ---
  H_obs[1:p, 1:p] <- crossprod(X, X * v)

  # --- Beta-theta and theta-theta blocks ---
  for (j in seq_len(J)) {
    idx_j <- which(group == j)
    v_j   <- v[idx_j]

    # Beta-theta_j cross block
    if (length(idx_j) == 1L) {
      bt_j <- v_j * X[idx_j, ]
    } else {
      bt_j <- colSums(X[idx_j, , drop = FALSE] * v_j)
    }
    H_obs[1:p, p + j] <- bt_j
    H_obs[p + j, 1:p] <- bt_j

    # Theta_j-theta_j diagonal
    H_obs[p + j, p + j] <- sum(v_j)
  }

  # --- Add prior precision to diagonal ---
  if (is.finite(beta_prior_prec) && beta_prior_prec > 0) {
    for (k in seq_len(p)) {
      H_obs[k, k] <- H_obs[k, k] + beta_prior_prec
    }
  }

  for (j in seq_len(J)) {
    H_obs[p + j, p + j] <- H_obs[p + j, p + j] + theta_prior_prec
  }

  H_obs
}


#' Safely invert a matrix with nearPD fallback
#'
#' Attempts \code{solve(H_obs)}. If that fails (singular or near-singular),
#' falls back to \code{Matrix::nearPD} to find the nearest positive-definite
#' matrix, then inverts that.
#'
#' @param H_obs Square numeric matrix.
#'
#' @return The inverse of H_obs (or its nearest PD approximation).
#'
#' @keywords internal
.safe_invert <- function(H_obs) {

  tryCatch(
    solve(H_obs),
    error = function(e) {
      if (requireNamespace("Matrix", quietly = TRUE)) {
        H_pd <- as.matrix(Matrix::nearPD(H_obs, keepDiag = TRUE)$mat)
        tryCatch(
          solve(H_pd),
          error = function(e2) {
            # If nearPD result is still singular, add a small ridge
            ridge <- max(abs(diag(H_pd))) * 1e-10
            solve(H_pd + ridge * diag(nrow(H_pd)))
          }
        )
      } else {
        stop("H_obs is singular and Matrix package is not available ",
             "for nearPD fallback.", call. = FALSE)
      }
    }
  )
}


#' Build the clustered score outer product matrix
#'
#' Constructs J_cluster = sum_g s_g s_g^T where s_g is the score vector
#' summed within PSU g.
#'
#' @param X Design matrix (N x p).
#' @param r Weighted residuals (length N).
#' @param psu Integer PSU indicator (1:G), length N.
#' @param group Integer group indicator (1:J), length N.
#' @param p Number of fixed-effect parameters.
#' @param J Number of groups.
#'
#' @return Symmetric (p+J) x (p+J) matrix.
#'
#' @keywords internal
.build_J_cluster <- function(X, r, psu, group, p, J) {

  d <- p + J
  G <- max(psu)
  J_cluster <- matrix(0, d, d)

  for (g in seq_len(G)) {
    idx_g <- which(psu == g)
    if (length(idx_g) == 0L) next

    s_g <- numeric(d)

    # Beta score component
    if (length(idx_g) == 1L) {
      s_g[1:p] <- r[idx_g] * X[idx_g, ]
    } else {
      s_g[1:p] <- colSums(X[idx_g, , drop = FALSE] * r[idx_g])
    }

    # Theta score component: accumulate within each group present in PSU g
    for (j_state in unique(group[idx_g])) {
      idx_gj <- idx_g[group[idx_g] == j_state]
      s_g[p + j_state] <- s_g[p + j_state] + sum(r[idx_gj])
    }

    J_cluster <- J_cluster + tcrossprod(s_g)
  }

  # Symmetrize
  J_cluster <- (J_cluster + t(J_cluster)) / 2

  J_cluster
}


#' Build the sandwich variance estimator
#'
#' Computes V_sand = H_obs_inv %*% J_cluster %*% H_obs_inv and symmetrizes.
#'
#' @param H_obs_inv Inverse of the observed information matrix (d x d).
#' @param J_cluster Clustered score outer product matrix (d x d).
#'
#' @return Symmetric (d x d) sandwich variance matrix.
#'
#' @keywords internal
.build_V_sand <- function(H_obs_inv, J_cluster) {

  V_sand <- H_obs_inv %*% J_cluster %*% H_obs_inv
  V_sand <- (V_sand + t(V_sand)) / 2

  V_sand
}
