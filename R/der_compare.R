###############################################################################
# der_compare.R
# Compare DER across clustering definitions
# ---------------------------------------------------------------------------
# Computes DER using different clustering (PSU) definitions and compares
# the results side by side.
###############################################################################

#' Compare DER Across Clustering Definitions
#'
#' Computes DER using different clustering (PSU) definitions
#' and compares the results. This is useful for comparing, e.g.,
#' state-level vs PSU-level clustering to assess sensitivity
#' of DER to the choice of primary sampling unit.
#'
#' @param x A draws matrix (S x d) or a \code{svyder} object.
#'   If a \code{svyder} object, the original draws and data
#'   are extracted automatically.
#' @param clusters A named list of PSU vectors to compare.
#'   Each element should be an integer vector of length N.
#' @param ... Additional arguments passed to [der_compute()].
#'   Required when \code{x} is a matrix: \code{y}, \code{X}, \code{group},
#'   \code{weights}, \code{family}, \code{sigma_theta}, etc.
#'
#' @return A \code{data.frame} with columns: \code{param}, \code{cluster_name},
#'   \code{der}.
#'
#' @seealso [der_compute()] for computing DER with a single clustering.
#' @family analysis
#'
#' @examples
#' data(nsece_demo)
#' # Compare DER using original PSU vs group-level clustering
#' comp <- der_compare(
#'   nsece_demo$draws,
#'   clusters = list(
#'     psu   = nsece_demo$psu,
#'     group = nsece_demo$group
#'   ),
#'   y = nsece_demo$y, X = nsece_demo$X,
#'   group = nsece_demo$group, weights = nsece_demo$weights,
#'   family = "binomial",
#'   sigma_theta = nsece_demo$sigma_theta,
#'   param_types = nsece_demo$param_types
#' )
#' head(comp)
#'
#' @export
der_compare <- function(x, clusters, ...) {

  stopifnot(is.list(clusters), length(clusters) >= 1L)

  if (is.null(names(clusters)) || any(names(clusters) == "")) {
    stop("'clusters' must be a named list. Each element names a clustering definition.",
         call. = FALSE)
  }

  # --- Extract draws and metadata from svyder object if needed ---
  if (is.svyder(x)) {
    draws_mat <- x$original_draws
    # User must supply the remaining arguments via ... since the svyder
    # object does not store y, X, group, weights separately.
    # We pass the draws matrix forward.
    x_input <- draws_mat
  } else {
    x_input <- x
  }

  cluster_names <- names(clusters)
  results_list  <- vector("list", length(cluster_names))

  for (k in seq_along(cluster_names)) {
    cname <- cluster_names[k]
    psu_k <- clusters[[k]]

    # Compute DER with this PSU definition
    svyder_k <- der_compute(x_input, ..., psu = psu_k)

    der_k       <- svyder_k$der
    param_names <- svyder_k$params

    results_list[[k]] <- data.frame(
      param        = param_names,
      cluster_name = rep(cname, length(param_names)),
      der          = as.numeric(der_k),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results_list)
}
