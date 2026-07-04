###############################################################################
# der_correct.R
# Selective correction of flagged parameters
# ---------------------------------------------------------------------------
# Two methods:
#   "block_cholesky" (default): joint Cholesky rescaling of the flagged
#       block, matching both marginal variances AND the covariance structure
#       of the sandwich target on that block. This is Algorithm 1 (CCC) of
#       the paper.
#   "marginal": per-parameter sqrt(DER) rescaling. Matches marginal
#       variances only; off-diagonal structure of the flagged block is NOT
#       matched. Coincides with block_cholesky when exactly one parameter
#       is flagged.
# Unflagged draws are left bitwise identical to the originals.
###############################################################################

#' Apply Selective Correction to Flagged Parameters
#'
#' Rescales the posterior draws of parameters flagged by [der_classify()]
#' so that their covariance matches the sandwich variance target, leaving
#' unflagged parameters untouched. The correction preserves posterior means.
#'
#' @details
#' With \code{method = "block_cholesky"} (the paper's Algorithm 1), the
#' flagged block F is transformed jointly: writing
#' \eqn{\Sigma_F = L_1 L_1^\top} and \eqn{V_F = L_2 L_2^\top} for the MCMC
#' and sandwich covariance of the flagged block, each draw is mapped through
#' \eqn{\phi^* = \hat\phi + L_2 L_1^{-1} (\phi - \hat\phi)}, so the corrected
#' block has covariance exactly \eqn{V_F}. With a single flagged parameter
#' this reduces to scaling by \eqn{\sqrt{V_{ii}/\Sigma_{ii}}}.
#'
#' With \code{method = "marginal"}, each flagged parameter is rescaled
#' independently by \eqn{\sqrt{V_{ii}/\Sigma_{ii}}}. This matches marginal
#' variances only and is provided for transparency and for reproducing
#' results computed this way; it is not a Cholesky correction when more
#' than one parameter is flagged.
#'
#' @param x A \code{svyder} object with classification (from
#'   [der_classify()]).
#' @param method \code{"block_cholesky"} (default) or \code{"marginal"}.
#'
#' @return A \code{svyder} object with \code{corrected_draws},
#'   \code{scale_factors} (marginal SD ratios of flagged parameters),
#'   and \code{correction_method} populated.
#'
#' @references
#' Lee, J., Williams, M. R., & Savitsky, T. D. (2026). Design Effect Ratios
#' for Bayesian Survey Models: A Diagnostic Framework for Identifying
#' Survey-Sensitive Parameters. \emph{Journal of Survey Statistics and
#' Methodology}. Submitted.
#'
#' Savitsky, T. D., & Williams, M. R. (2022). Pseudo Bayesian mixed models
#' under informative sampling. \emph{Journal of Official Statistics}, 38(4),
#' 901--928. \doi{10.2478/jos-2022-0039}
#'
#' Williams, M. R., & Savitsky, T. D. (2021). Uncertainty estimation for
#' pseudo-Bayesian inference under complex sampling. \emph{International
#' Statistical Review}, 89(1), 72--107. \doi{10.1111/insr.12376}
#'
#' Higham, N. J. (2002). Computing the nearest correlation matrix---a problem
#' from finance. \emph{IMA Journal of Numerical Analysis}, 22(3), 329--343.
#' \doi{10.1093/imanum/22.3.329}
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
#'   cluster = nsece_demo$psu, family = "binomial",
#'   sigma_theta = nsece_demo$sigma_theta,
#'   param_types = nsece_demo$param_types
#' )
#' result <- der_classify(result, tau = 1.2, verbose = FALSE)
#' result <- der_correct(result)
#'
#' @export
der_correct <- function(x, method = c("block_cholesky", "marginal")) {

  if (!is.svyder(x)) {
    stop("'x' must be a svyder object.", call. = FALSE)
  }

  if (!("flagged" %in% names(x$classification))) {
    stop("Classification not found. Run der_classify() first.", call. = FALSE)
  }

  if (identical(method, "cholesky")) {
    stop("method = 'cholesky' was ambiguous and has been retired.\n",
         "  * 'block_cholesky': joint correction of the flagged block ",
         "(the paper's Algorithm 1; new default)\n",
         "  * 'marginal': the previous per-parameter sqrt(DER) rescaling",
         call. = FALSE)
  }
  method <- match.arg(method)

  d          <- length(x$der)
  draws_all  <- x$original_draws
  diag_V     <- diag(x$V_sand)
  diag_mcmc  <- diag(x$sigma_mcmc)
  flagged    <- x$classification$flagged
  point_est  <- colMeans(draws_all)

  scale_factors   <- rep(1.0, d)
  draws_corrected <- draws_all

  F_idx     <- which(flagged)
  n_flagged <- length(F_idx)

  if (n_flagged > 0) {
    # Marginal SD ratios (reported for both methods; the block method also
    # matches these marginals exactly)
    scale_factors[F_idx] <- sqrt(diag_V[F_idx] / diag_mcmc[F_idx])

    if (method == "marginal" || n_flagged == 1L) {
      for (i in F_idx) {
        draws_corrected[, i] <- point_est[i] +
          scale_factors[i] * (draws_all[, i] - point_est[i])
      }
    } else {
      Sigma_F <- x$sigma_mcmc[F_idx, F_idx, drop = FALSE]
      V_F     <- x$V_sand[F_idx, F_idx, drop = FALSE]

      L1 <- tryCatch(
        t(chol(Sigma_F)),
        error = function(e) {
          stop("MCMC covariance of the flagged block is not positive ",
               "definite; cannot apply block correction.", call. = FALSE)
        }
      )
      L2 <- tryCatch(
        t(chol(V_F)),
        error = function(e) {
          if (!requireNamespace("Matrix", quietly = TRUE)) {
            stop("Sandwich covariance of the flagged block is not positive ",
                 "definite and the Matrix package is unavailable for a ",
                 "nearPD fallback.", call. = FALSE)
          }
          warning("Sandwich covariance of the flagged block is not positive ",
                  "definite (this can occur when the number of clusters is ",
                  "smaller than the flagged block); using the nearest ",
                  "positive-definite matrix.", call. = FALSE)
          t(chol(as.matrix(Matrix::nearPD(V_F)$mat)))
        }
      )

      A <- L2 %*% solve(L1)
      centered <- sweep(draws_all[, F_idx, drop = FALSE], 2L, point_est[F_idx])
      draws_corrected[, F_idx] <- sweep(centered %*% t(A), 2L,
                                        point_est[F_idx], FUN = "+")
    }
  }

  x$corrected_draws   <- draws_corrected
  x$scale_factors     <- scale_factors
  x$correction_method <- method

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
#' @references
#' Lee, J., Williams, M. R., & Savitsky, T. D. (2026). Design Effect Ratios
#' for Bayesian Survey Models: A Diagnostic Framework for Identifying
#' Survey-Sensitive Parameters. \emph{Journal of Survey Statistics and
#' Methodology}. Submitted.
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
#'   cluster = nsece_demo$psu, family = "binomial",
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
