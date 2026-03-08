###############################################################################
# extract-posterior.R
# posterior package integration for the svyder package
# ---------------------------------------------------------------------------
# Provides extract_draws methods for the posterior package's draws formats:
#   - draws_matrix
#   - draws_df
#   - draws_array
#   - draws_list
#   - draws_rvars
#
# These are conversion methods that allow users to pass posterior draws
# objects directly to svyder functions.
###############################################################################

#' @rdname extract_draws
#' @export
extract_draws.draws_matrix <- function(x, ...) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for extract_draws.draws_matrix(). ",
         "Install it with: install.packages('posterior')",
         call. = FALSE)
  }
  draws_mat <- as.matrix(x)
  # Remove .chain, .iteration, .draw columns if present
  meta_cols <- c(".chain", ".iteration", ".draw")
  keep <- !(colnames(draws_mat) %in% meta_cols)
  draws_mat <- draws_mat[, keep, drop = FALSE]

  list(draws = draws_mat)
}

#' @rdname extract_draws
#' @export
extract_draws.draws_df <- function(x, ...) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for extract_draws.draws_df(). ",
         "Install it with: install.packages('posterior')",
         call. = FALSE)
  }
  # Convert to draws_matrix first, then to plain matrix
  draws_mat <- as.matrix(posterior::as_draws_matrix(x))
  # Remove metadata columns if present
  meta_cols <- c(".chain", ".iteration", ".draw")
  keep <- !(colnames(draws_mat) %in% meta_cols)
  draws_mat <- draws_mat[, keep, drop = FALSE]

  list(draws = draws_mat)
}

#' @rdname extract_draws
#' @export
extract_draws.draws_array <- function(x, ...) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for extract_draws.draws_array(). ",
         "Install it with: install.packages('posterior')",
         call. = FALSE)
  }
  # Convert to draws_matrix (merging chains), then to plain matrix
  draws_mat <- as.matrix(posterior::as_draws_matrix(x))
  meta_cols <- c(".chain", ".iteration", ".draw")
  keep <- !(colnames(draws_mat) %in% meta_cols)
  draws_mat <- draws_mat[, keep, drop = FALSE]

  list(draws = draws_mat)
}

#' @rdname extract_draws
#' @export
extract_draws.draws_list <- function(x, ...) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for extract_draws.draws_list(). ",
         "Install it with: install.packages('posterior')",
         call. = FALSE)
  }
  draws_mat <- as.matrix(posterior::as_draws_matrix(x))
  meta_cols <- c(".chain", ".iteration", ".draw")
  keep <- !(colnames(draws_mat) %in% meta_cols)
  draws_mat <- draws_mat[, keep, drop = FALSE]

  list(draws = draws_mat)
}

#' @rdname extract_draws
#' @export
extract_draws.draws_rvars <- function(x, ...) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for extract_draws.draws_rvars(). ",
         "Install it with: install.packages('posterior')",
         call. = FALSE)
  }
  draws_mat <- as.matrix(posterior::as_draws_matrix(x))
  meta_cols <- c(".chain", ".iteration", ".draw")
  keep <- !(colnames(draws_mat) %in% meta_cols)
  draws_mat <- draws_mat[, keep, drop = FALSE]

  list(draws = draws_mat)
}
