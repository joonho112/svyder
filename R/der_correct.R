###############################################################################
# der_correct.R
# Selective Cholesky correction for flagged parameters
# ---------------------------------------------------------------------------
# Rescales posterior draws for flagged parameters so that their marginal
# variance matches the sandwich variance. Unflagged draws are left bitwise
# identical to the originals.
###############################################################################

#' Apply Selective Cholesky Correction
#'
#' For each parameter flagged by [der_classify()], rescales the
#' posterior draws so that marginal variance matches the sandwich variance
#' estimate. The correction preserves the posterior mean: draws are centered,
#' scaled by \code{sqrt(V_sand[i,i] / sigma_mcmc[i,i])}, then re-centered.
#'
#' Unflagged parameters retain their original draws without any modification.
#'
#' @param x A \code{svyder} object with classification (from
#'   [der_classify()]).
#' @param method Correction method (default \code{"cholesky"}). Currently only
#'   \code{"cholesky"} is supported.
#'
#' @return A \code{svyder} object with \code{corrected_draws},
#'   \code{scale_factors}, and \code{original_draws} populated.
#'
#' @seealso [der_classify()] for flagging parameters, [der_compute()] for
#'   computing DER, [as.matrix.svyder()] for extracting corrected draws.
#' @family core-pipeline
#'
#' @examples
#' data(nsece_demo)
#' result <- der_compute(
#'   nsece_demo$draws,
#'   y = nsece_demo$y, X = nsece_demo$X,
#'   group = nsece_demo$group, weights = nsece_demo$weights,
#'   psu = nsece_demo$psu, family = "binomial",
#'   sigma_theta = nsece_demo$sigma_theta,
#'   param_types = nsece_demo$param_types
#' )
#' result <- der_classify(result, tau = 1.2, verbose = FALSE)
#' result <- der_correct(result)
#'
#' @export
der_correct <- function(x, method = "cholesky") {

  if (!is.svyder(x)) {
    stop("'x' must be a svyder object.", call. = FALSE)
  }

  if (!("flagged" %in% names(x$classification))) {
    stop("Classification not found. Run der_classify() first.", call. = FALSE)
  }

  match.arg(method, choices = "cholesky")

  d            <- length(x$der)
  draws_all    <- x$original_draws
  diag_V       <- diag(x$V_sand)
  diag_mcmc    <- diag(x$sigma_mcmc)
  flagged      <- x$classification$flagged
  point_est    <- colMeans(draws_all)

  scale_factors     <- rep(1.0, d)
  draws_corrected   <- draws_all

  n_flagged <- sum(flagged)

  if (n_flagged > 0) {
    for (i in seq_len(d)) {
      if (flagged[i]) {
        sf <- sqrt(diag_V[i] / diag_mcmc[i])
        scale_factors[i] <- sf
        draws_corrected[, i] <- point_est[i] +
          sf * (draws_all[, i] - point_est[i])
      }
    }
  }

  x$corrected_draws <- draws_corrected
  x$scale_factors   <- scale_factors

  x
}


#' Extract Draws Matrix from a svyder Object
#'
#' Returns the corrected draws if available, otherwise the original draws.
#' This method allows \code{svyder} objects to be used wherever a numeric
#' matrix of posterior draws is expected.
#'
#' @param x A \code{svyder} object.
#' @param ... Ignored.
#'
#' @return A numeric matrix of posterior draws (S x d).
#'
#' @seealso [der_correct()] for applying the correction.
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
#' draws <- as.matrix(result)
#' dim(draws)
#'
#' @export
as.matrix.svyder <- function(x, ...) {
  if (!is.null(x$corrected_draws)) x$corrected_draws else x$original_draws
}
