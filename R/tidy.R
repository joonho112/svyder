###############################################################################
# tidy.R
# Tidy and glance methods for svyder objects
# ---------------------------------------------------------------------------
# Returns tidy data frames following broom conventions.
###############################################################################

#' Tidy a svyder Object
#'
#' Returns a one-row-per-parameter data frame with DER diagnostics
#' in a tidy format compatible with broom conventions. Includes
#' posterior summaries, classification tier, and correction scale factor.
#'
#' @param x A \code{svyder} object.
#' @param ... Additional arguments (unused).
#'
#' @return A \code{data.frame} with one row per parameter and columns:
#'   \describe{
#'     \item{term}{Parameter name.}
#'     \item{estimate}{Posterior mean (from original draws).}
#'     \item{std.error}{Posterior standard deviation (from original draws).}
#'     \item{der}{Design Effect Ratio.}
#'     \item{tier}{Three-tier classification (if classified).}
#'     \item{action}{Action label: \code{"CORRECT"} or \code{"retain"}
#'       (if classified).}
#'     \item{flagged}{Logical; whether the parameter is flagged for
#'       correction (if classified).}
#'     \item{scale_factor}{Cholesky scale factor applied to this parameter.}
#'   }
#'
#' @seealso [glance.svyder()] for model-level summaries,
#'   [print.svyder()] for console output.
#' @family svyder-methods
#'
#' @examples
#' data(nsece_demo)
#' result <- der_diagnose(
#'   nsece_demo$draws,
#'   y = nsece_demo$y, X = nsece_demo$X,
#'   group = nsece_demo$group, weights = nsece_demo$weights,
#'   psu = nsece_demo$psu, family = "binomial",
#'   sigma_theta = nsece_demo$sigma_theta,
#'   param_types = nsece_demo$param_types
#' )
#' head(tidy.svyder(result))
#'
#' @export
tidy.svyder <- function(x, ...) {

  d <- length(x$der)
  has_class <- "flagged" %in% names(x$classification)

  # Posterior summaries from original draws
  if (!is.null(x$original_draws)) {
    estimate  <- colMeans(x$original_draws)
    std_error <- apply(x$original_draws, 2, stats::sd)
  } else {
    estimate  <- rep(NA_real_, d)
    std_error <- rep(NA_real_, d)
  }

  out <- data.frame(
    term         = x$params,
    estimate     = estimate,
    std.error    = std_error,
    der          = as.numeric(x$der),
    stringsAsFactors = FALSE
  )

  if (has_class) {
    out$tier        <- x$classification$tier
    out$action      <- x$classification$action
    out$flagged     <- x$classification$flagged
  } else {
    out$tier        <- NA_character_
    out$action      <- NA_character_
    out$flagged     <- NA
  }

  out$scale_factor <- x$scale_factors

  out
}


#' Glance at a svyder Object
#'
#' Returns a one-row data frame with model-level summary statistics
#' following broom conventions. Provides an overview of the DER
#' diagnostic results.
#'
#' @param x A \code{svyder} object.
#' @param ... Additional arguments (unused).
#'
#' @return A \code{data.frame} with one row and columns:
#'   \describe{
#'     \item{n_params}{Total number of parameters.}
#'     \item{n_flagged}{Number of parameters flagged for correction.}
#'     \item{pct_flagged}{Percentage of parameters flagged.}
#'     \item{tau}{Classification threshold.}
#'     \item{family}{Model family (\code{"binomial"} or \code{"gaussian"}).}
#'     \item{n_obs}{Number of observations.}
#'     \item{n_groups}{Number of groups/clusters.}
#'     \item{mean_deff}{Mean per-group design effect.}
#'     \item{mean_B}{Mean per-group shrinkage factor.}
#'     \item{der_min}{Minimum DER value.}
#'     \item{der_max}{Maximum DER value.}
#'   }
#'
#' @seealso [tidy.svyder()] for per-parameter summaries.
#' @family svyder-methods
#'
#' @examples
#' data(nsece_demo)
#' result <- der_diagnose(
#'   nsece_demo$draws,
#'   y = nsece_demo$y, X = nsece_demo$X,
#'   group = nsece_demo$group, weights = nsece_demo$weights,
#'   psu = nsece_demo$psu, family = "binomial",
#'   sigma_theta = nsece_demo$sigma_theta,
#'   param_types = nsece_demo$param_types
#' )
#' glance.svyder(result)
#'
#' @export
glance.svyder <- function(x, ...) {

  d <- length(x$der)
  has_class <- "flagged" %in% names(x$classification)

  n_flagged <- if (has_class) sum(x$classification$flagged) else NA_integer_
  pct_flagged <- if (has_class) n_flagged / d * 100 else NA_real_

  data.frame(
    n_params    = d,
    n_flagged   = n_flagged,
    pct_flagged = pct_flagged,
    tau         = x$tau,
    family      = x$family,
    n_obs       = x$n_obs,
    n_groups    = x$n_groups,
    mean_deff   = mean(x$deff_j),
    mean_B      = mean(x$B_j),
    der_min     = min(x$der),
    der_max     = max(x$der),
    stringsAsFactors = FALSE
  )
}
