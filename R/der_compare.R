###############################################################################
# der_compare.R
# Compare DER across variance targets (aggregation units)
# ---------------------------------------------------------------------------
# The sandwich aggregation unit is part of the variance-target definition;
# this function reports DER under several declared targets side by side.
###############################################################################

#' Compare DER Across Variance Targets
#'
#' Computes DER under different sandwich aggregation units (e.g., the
#' design PSU vs the model group) and compares the results side by side.
#' Because the aggregation unit is part of the variance-target definition,
#' reporting both is the recommended practice whenever the design PSUs and
#' the model groups differ.
#'
#' When \code{x} is a \code{svyder} object from svyder >= 0.2.0, the stored
#' data slots are reused: only the meat is rebuilt per target (the bread and
#' the MCMC covariance are target-invariant), so the comparison is fast and
#' requires no re-supply of \code{y}/\code{X}/\code{weights}.
#' If the stored object used \code{cluster_strata} to declare selected-but-empty
#' clusters, that declared universe is preserved when a comparison target is
#' exactly the original stored target.
#'
#' @param x A \code{svyder} object (preferred) or a draws matrix (S x d).
#' @param clusters A named list of cluster-id vectors (each length N), one
#'   per target to compare.
#' @param ... For the matrix path only: arguments passed to [der_compute()]
#'   (\code{y}, \code{X}, \code{group}, \code{weights}, \code{family},
#'   \code{sigma_theta}, ...).
#'
#' @return A \code{data.frame} with columns: \code{param},
#'   \code{cluster_name}, \code{der}.
#'
#' @references
#' Lee, J., Williams, M. R., & Savitsky, T. D. (2026). Design Effect Ratios
#' for Bayesian Survey Models: A Diagnostic Framework for Identifying
#' Survey-Sensitive Parameters. \emph{Journal of Survey Statistics and
#' Methodology}. Submitted.
#'
#' @seealso [der_compute()] for computing DER under a single target.
#' @family analysis
#'
#' @examples
#' data(nsece_demo)
#' result <- der_compute(
#'   nsece_demo$draws,
#'   y = nsece_demo$y, X = nsece_demo$X,
#'   group = nsece_demo$group, weights = nsece_demo$weights,
#'   cluster = nsece_demo$psu, family = "binomial",
#'   sigma_theta = nsece_demo$sigma_theta,
#'   param_types = nsece_demo$param_types
#' )
#' comp <- der_compare(result, clusters = list(
#'   design_psu  = nsece_demo$psu,
#'   model_group = nsece_demo$group
#' ))
#' head(comp)
#'
#' @export
der_compare <- function(x, clusters, ...) {

  stopifnot(is.list(clusters), length(clusters) >= 1L)

  if (is.null(names(clusters)) || any(names(clusters) == "")) {
    stop("'clusters' must be a named list. Each element names a variance target.",
         call. = FALSE)
  }

  cluster_names <- names(clusters)
  results_list  <- vector("list", length(cluster_names))

  if (is.svyder(x) && !is.null(x$data)) {
    # --- Fast path: rebuild the meat only ---
    H_obs_inv <- .safe_invert(x$H_obs)
    diag_mcmc <- diag(x$sigma_mcmc)
    J <- x$n_groups
    p <- length(x$der) - J

    for (k in seq_along(cluster_names)) {
      cname     <- cluster_names[k]
      cluster_k <- clusters[[k]]

      if (length(cluster_k) != nrow(x$data$X)) {
        stop(sprintf(
          "clusters[['%s']] has length %d but the data have %d observations.",
          cname, length(cluster_k), nrow(x$data$X)
        ), call. = FALSE)
      }
      raw_cluster_int <- suppressWarnings({
        if (is.factor(cluster_k)) {
          as.integer(as.character(cluster_k))
        } else {
          as.integer(cluster_k)
        }
      })
      raw_cluster_ok <- !anyNA(raw_cluster_int)
      if (raw_cluster_ok && is.numeric(cluster_k)) {
        raw_cluster_ok <- all(as.numeric(cluster_k) == raw_cluster_int)
      }

      same_declared_target <- !is.null(x$data$cluster_strata) &&
        raw_cluster_ok &&
        length(raw_cluster_int) == length(x$data$cluster) &&
        identical(as.integer(raw_cluster_int), as.integer(x$data$cluster))

      if (same_declared_target) {
        cluster_int <- x$data$cluster
        cluster_strata_k <- x$data$cluster_strata
      } else {
        cluster_int <- as.integer(as.factor(cluster_k))
        cluster_strata_k <- NULL
      }

      build_stratified <- function() {
        .build_J_cluster(x$data$X, x$data$r, cluster_int, x$data$group,
                         p, J,
                         strata         = x$data$strata,
                         cluster_strata = cluster_strata_k,
                         center         = x$target$center,
                         df_correct     = x$target$df_correct)
      }

      # Stored strata are used when the new clusters nest within them;
      # otherwise the target is computed unstratified with a message.
      J_k <- if (same_declared_target) {
        build_stratified()
      } else {
        tryCatch(
          build_stratified(),
          error = function(e) {
            message(sprintf(
              "Target '%s': clusters do not nest within the stored strata; computing unstratified.",
              cname
            ))
            .build_J_cluster(x$data$X, x$data$r, cluster_int, x$data$group,
                             p, J,
                             strata         = NULL,
                             cluster_strata = NULL,
                             center         = x$target$center,
                             df_correct     = x$target$df_correct)
          }
        )
      }

      V_k   <- .build_V_sand(H_obs_inv, J_k)
      der_k <- diag(V_k) / diag_mcmc

      results_list[[k]] <- data.frame(
        param        = x$params,
        cluster_name = rep(cname, length(x$params)),
        der          = as.numeric(der_k),
        stringsAsFactors = FALSE
      )
    }

    return(do.call(rbind, results_list))
  }

  # --- Matrix path: full recomputation per target ---
  x_input <- if (is.svyder(x)) x$original_draws else x

  for (k in seq_along(cluster_names)) {
    cname     <- cluster_names[k]
    cluster_k <- clusters[[k]]

    svyder_k <- der_compute(x_input, ..., cluster = cluster_k)

    results_list[[k]] <- data.frame(
      param        = svyder_k$params,
      cluster_name = rep(cname, length(svyder_k$params)),
      der          = as.numeric(svyder_k$der),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results_list)
}
