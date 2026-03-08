###############################################################################
# utils.R
# Input validation and helper utilities for the svyder package
# ---------------------------------------------------------------------------
# All functions are internal (not exported).
###############################################################################

# Comprehensive input validation with informative error messages
#
# @param y Numeric/integer vector of responses (length N)
# @param X Numeric matrix of predictors (N x p)
# @param group Integer vector of group indicators (length N)
# @param weights Numeric vector of survey weights (length N)
# @param family Character string: "binomial" or "gaussian"
# @param draws_beta Numeric matrix of posterior draws for beta (M x p)
# @param draws_theta Numeric matrix of posterior draws for theta (M x J)
# @return Invisible TRUE if all checks pass; otherwise stops with error
.validate_inputs <- function(y, X, group, weights, family, draws_beta,
                             draws_theta) {
  .stop_msg <- function(msg) {
    if (requireNamespace("cli", quietly = TRUE)) {
      cli::cli_abort(msg)
    } else {
      stop(msg, call. = FALSE)
    }
  }

  # --- y ---
  if (!is.numeric(y) && !is.integer(y)) {
    .stop_msg("'y' must be a numeric or integer vector.")
  }
  N <- length(y)
  if (N == 0L) {
    .stop_msg("'y' has length 0.")
  }

  # --- X ---
  if (!is.matrix(X)) {
    .stop_msg("'X' must be a matrix.")
  }
  if (nrow(X) != N) {
    .stop_msg(sprintf(
      "Dimension mismatch: length(y) = %d but nrow(X) = %d.",
      N, nrow(X)
    ))
  }
  p <- ncol(X)

  # --- group ---
  if (!is.numeric(group) && !is.integer(group)) {
    .stop_msg("'group' must be a numeric or integer vector.")
  }
  if (length(group) != N) {
    .stop_msg(sprintf(
      "Dimension mismatch: length(y) = %d but length(group) = %d.",
      N, length(group)
    ))
  }
  group_int <- as.integer(group)
  J <- max(group_int)
  expected_groups <- seq_len(J)
  actual_groups <- sort(unique(group_int))
  if (!identical(actual_groups, expected_groups)) {
    .stop_msg(sprintf(
      "'group' must contain sequential integers from 1 to J (max = %d). Missing groups: %s.",
      J,
      paste(setdiff(expected_groups, actual_groups), collapse = ", ")
    ))
  }

  # --- weights ---
  if (!is.numeric(weights)) {
    .stop_msg("'weights' must be a numeric vector.")
  }
  if (length(weights) != N) {
    .stop_msg(sprintf(
      "Dimension mismatch: length(y) = %d but length(weights) = %d.",
      N, length(weights)
    ))
  }
  if (any(weights <= 0)) {
    n_bad <- sum(weights <= 0)
    .stop_msg(sprintf(
      "All survey weights must be positive. Found %d non-positive weight(s).",
      n_bad
    ))
  }

  # --- family ---
  valid_families <- c("binomial", "gaussian")
  if (!is.character(family) || length(family) != 1L || !(family %in% valid_families)) {
    .stop_msg(sprintf(
      "'family' must be one of: %s. Got: '%s'.",
      paste(sprintf("'%s'", valid_families), collapse = ", "),
      as.character(family)
    ))
  }

  # --- draws_beta ---
  if (!is.matrix(draws_beta)) {
    .stop_msg("'draws_beta' must be a matrix.")
  }
  if (ncol(draws_beta) != p) {
    .stop_msg(sprintf(
      "Dimension mismatch: ncol(X) = %d but ncol(draws_beta) = %d.",
      p, ncol(draws_beta)
    ))
  }
  M <- nrow(draws_beta)

  # --- draws_theta ---
  if (!is.matrix(draws_theta)) {
    .stop_msg("'draws_theta' must be a matrix.")
  }
  if (ncol(draws_theta) != J) {
    .stop_msg(sprintf(
      "Dimension mismatch: max(group) = %d but ncol(draws_theta) = %d.",
      J, ncol(draws_theta)
    ))
  }
  if (nrow(draws_theta) != M) {
    .stop_msg(sprintf(
      "Dimension mismatch: nrow(draws_beta) = %d but nrow(draws_theta) = %d.",
      M, nrow(draws_theta)
    ))
  }

  invisible(TRUE)
}


# Validate svyder object structure
#
# @param x An object to check
# @return Invisible TRUE if valid; otherwise stops with error
.validate_svyder <- function(x) {
  .stop_msg <- function(msg) {
    if (requireNamespace("cli", quietly = TRUE)) {
      cli::cli_abort(msg)
    } else {
      stop(msg, call. = FALSE)
    }
  }

  if (!inherits(x, "svyder")) {
    .stop_msg("Object must be of class 'svyder'.")
  }

  required_components <- c("der", "classification", "V_sand", "sigma_mcmc",
                           "diagnostics")
  missing <- setdiff(required_components, names(x))
  if (length(missing) > 0L) {
    .stop_msg(sprintf(
      "svyder object is missing required component(s): %s.",
      paste(sprintf("'%s'", missing), collapse = ", ")
    ))
  }

  invisible(TRUE)
}


# Compute credible intervals from posterior draws
#
# @param draws Numeric matrix of posterior draws (M x d)
# @param prob Numeric scalar, credible interval probability (default 0.95)
# @return Numeric matrix with columns "lower" and "upper" (d x 2)
.compute_ci <- function(draws, prob = 0.95) {
  if (!is.matrix(draws)) {
    draws <- matrix(draws, ncol = 1L)
  }

  alpha <- (1 - prob) / 2
  probs <- c(alpha, 1 - alpha)

  ci <- t(apply(draws, 2L, quantile, probs = probs))
  colnames(ci) <- c("lower", "upper")
  ci
}


# Format numeric values as percentage strings
#
# @param x Numeric vector
# @param digits Integer, number of decimal places (default 1)
# @return Character vector of formatted percentage strings
.format_pct <- function(x, digits = 1L) {
  sprintf(paste0("%.", digits, "f%%"), x * 100)
}
