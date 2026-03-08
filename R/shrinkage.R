###############################################################################
# shrinkage.R
# Per-group design effect and shrinkage computations for the svyder package
# ---------------------------------------------------------------------------
# All functions are internal (not exported).
###############################################################################

# Per-group Kish design effect (DEFF)
# Ported from standalone compute_der.R lines 173-180
#
# DEFF_j = n_j * sum(w_j^2) / (sum(w_j))^2
#
# @param group Integer vector of group indicators (1:J), length N
# @param w Numeric vector of survey weights, length N
# @return Named numeric vector of per-group DEFF values (length J)
.compute_deff_j <- function(group, w) {
  J <- max(group)
  deff_j <- numeric(J)

  for (j in seq_len(J)) {
    idx_j <- which(group == j)
    w_j <- w[idx_j]
    n_j <- length(idx_j)
    deff_j[j] <- n_j * sum(w_j^2) / (sum(w_j))^2
  }

  names(deff_j) <- paste0("group_", seq_len(J))
  deff_j
}

# Per-group shrinkage factor B_j
# Ported from standalone compute_der.R lines 181-183
#
# B_j = sigma2 / (sigma2 + V_tilde_j)
# where V_tilde_j = 1 / sum(w_j * wt_j)
# and wt = working weights WITHOUT survey weights (e.g., mu*(1-mu) for binomial)
#
# @param group Integer vector of group indicators (1:J), length N
# @param w Numeric vector of survey weights, length N
# @param wt Numeric vector of working weights WITHOUT survey weights, length N
#   For binomial: wt_i = mu_i * (1 - mu_i)
#   For gaussian: wt_i = 1 / sigma_e^2
# @param sigma2 Numeric scalar, random effect variance (sigma_theta^2)
# @return Named numeric vector of per-group shrinkage factors (length J)
.compute_B_j <- function(group, w, wt, sigma2) {
  J <- max(group)
  B_j <- numeric(J)

  for (j in seq_len(J)) {
    idx_j <- which(group == j)
    info_j <- sum(w[idx_j] * wt[idx_j])
    V_tilde_j <- 1 / info_j
    B_j[j] <- sigma2 / (sigma2 + V_tilde_j)
  }

  names(B_j) <- paste0("group_", seq_len(J))
  B_j
}

# Finite-J correction factor (kappa)
#
# kappa_j = (J - 1) * (1 - B_j) / (J * (1 - B_j) + B_j)
#
# When B = 1 (complete shrinkage): kappa = 0
# When B -> 0 (no shrinkage): kappa -> (J-1)/J
#
# @param B_j Numeric vector of shrinkage factors
# @param J Integer, total number of groups
# @return Numeric vector of correction factors (same length as B_j)
.compute_kappa_j <- function(B_j, J) {
  (J - 1) * (1 - B_j) / (J * (1 - B_j) + B_j)
}
