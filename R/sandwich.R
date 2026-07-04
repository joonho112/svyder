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


#' Build the clustered score outer product matrix (meat)
#'
#' Constructs the meat matrix from cluster-level score totals. The general
#' (default) form is the stratified, centered, DF-corrected estimator
#' \deqn{J = \sum_h \frac{C_h}{C_h - 1} \sum_c (t_{hc} - \bar t_h)(t_{hc} - \bar t_h)^\top,}
#' where \eqn{t_{hc}} is the score total of cluster \eqn{c} in stratum
#' \eqn{h}, \eqn{C_h} the number of clusters in stratum \eqn{h}, and
#' \eqn{\bar t_h} the stratum mean of cluster totals. Centering matters at
#' the posterior-mean plug-in point because the likelihood score totals do
#' not sum to zero there (the prior gradient does not vanish).
#'
#' With \code{strata = NULL}, \code{center = FALSE}, \code{df_correct = FALSE}
#' this reduces to the uncentered single-stratum form
#' \eqn{J = \sum_c t_c t_c^\top} used by the v1 pipeline (kept for exact
#' reproduction of previously published results).
#'
#' Strata containing a single cluster cannot be centered; their uncentered
#' contribution is used with a warning.
#'
#' @param X Design matrix (N x p).
#' @param r Weighted score residuals (length N).
#' @param cluster Integer cluster indicator (1:G), length N. This is the
#'   sandwich aggregation unit (design PSU or model group).
#' @param group Integer group indicator (1:J), length N.
#' @param p Number of fixed-effect parameters.
#' @param J Number of groups.
#' @param strata Optional integer stratum indicator (1:H), length N.
#'   \code{NULL} treats the design as a single stratum.
#' @param cluster_strata Optional integer vector giving the stratum of
#'   every cluster in the declared design universe, indexed by cluster id
#'   (length G, the universe size). Supplying it declares the cluster
#'   universe explicitly: cluster ids in \code{1:G} that appear in no
#'   observation -- e.g., selected PSUs whose second-stage Poisson sample
#'   is empty -- remain in the meat as zero score-total clusters, entering
#'   the stratum centering and the C_h/(C_h-1) correction. Without it, the
#'   universe is inferred from the observed ids, which drops empty
#'   clusters. Requires \code{strata}.
#' @param center Center cluster score totals within strata (default TRUE).
#' @param df_correct Apply the C_h / (C_h - 1) stratum correction
#'   (default TRUE).
#'
#' @return Symmetric (p+J) x (p+J) matrix.
#'
#' @keywords internal
.build_J_cluster <- function(X, r, cluster, group, p, J,
                             strata = NULL, cluster_strata = NULL,
                             center = TRUE, df_correct = TRUE) {

  d <- p + J
  G <- if (!is.null(cluster_strata)) length(cluster_strata) else max(cluster)

  # --- Cluster-to-stratum map (each cluster must lie in one stratum) ---
  if (!is.null(cluster_strata)) {
    if (is.null(strata)) {
      stop("'cluster_strata' requires 'strata'.", call. = FALSE)
    }
    if (any(cluster < 1L | cluster > G)) {
      stop("With 'cluster_strata', cluster ids must lie in 1:length(cluster_strata).",
           call. = FALSE)
    }
    strata_of_cluster <- as.integer(cluster_strata)
    bad <- strata != strata_of_cluster[cluster]
    if (any(bad)) {
      stop(sprintf(
        "'strata' disagrees with 'cluster_strata' at %d observation(s).",
        sum(bad)), call. = FALSE)
    }
  } else if (is.null(strata)) {
    strata_of_cluster <- rep(1L, G)
  } else {
    tab <- unique(data.frame(cluster = cluster, strata = strata))
    if (nrow(tab) != G) {
      stop("Each cluster must belong to exactly one stratum.", call. = FALSE)
    }
    strata_of_cluster <- integer(G)
    strata_of_cluster[tab$cluster] <- tab$strata
  }

  # --- Cluster-level score totals: T is G x d ---
  T_mat <- matrix(0, G, d)
  for (g in seq_len(G)) {
    idx_g <- which(cluster == g)
    if (length(idx_g) == 0L) next

    s_g <- numeric(d)

    # Beta score component
    if (length(idx_g) == 1L) {
      s_g[1:p] <- r[idx_g] * X[idx_g, ]
    } else {
      s_g[1:p] <- colSums(X[idx_g, , drop = FALSE] * r[idx_g])
    }

    # Theta score component: accumulate within each group present in cluster g
    for (j_grp in unique(group[idx_g])) {
      idx_gj <- idx_g[group[idx_g] == j_grp]
      s_g[p + j_grp] <- s_g[p + j_grp] + sum(r[idx_gj])
    }

    T_mat[g, ] <- s_g
  }

  # --- Stratum-wise (centered, DF-corrected) accumulation ---
  J_cluster <- matrix(0, d, d)
  singleton_strata <- 0L

  for (h in unique(strata_of_cluster)) {
    idx_h <- which(strata_of_cluster == h)
    C_h   <- length(idx_h)
    T_h   <- T_mat[idx_h, , drop = FALSE]

    if (center && C_h > 1L) {
      T_h <- sweep(T_h, 2L, colMeans(T_h))
    } else if (center && C_h == 1L) {
      singleton_strata <- singleton_strata + 1L
    }

    mult <- if (df_correct && C_h > 1L) C_h / (C_h - 1) else 1
    J_cluster <- J_cluster + mult * crossprod(T_h)
  }

  if (singleton_strata > 0L) {
    warning(sprintf(
      "%d stratum/strata contain a single cluster; their contribution is uncentered.",
      singleton_strata
    ), call. = FALSE)
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
