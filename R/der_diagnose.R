###############################################################################
# der_diagnose.R
# All-in-one DER diagnostic pipeline
# ---------------------------------------------------------------------------
# Convenience wrapper that runs compute -> classify -> correct in sequence.
###############################################################################

#' All-in-One DER Diagnostic
#'
#' Runs the full DER pipeline in a single call: [der_compute()] to obtain
#' DER values, [der_classify()] to assign tiers and flag parameters, and
#' optionally [der_correct()] to apply selective correction. This is the
#' recommended entry point for most users.
#'
#' @inheritParams der_compute
#' @param tau Classification threshold (default 1.2).
#' @param correct Apply correction to flagged parameters (default \code{TRUE}).
#'
#' @return A fully processed \code{svyder} object with DER values,
#'   classification, and (optionally) corrected draws.
#'
#' @seealso [der_compute()], [der_classify()], [der_correct()] for the
#'   individual pipeline steps. [tidy.svyder()] and [glance.svyder()] for
#'   tidy summaries.
#' @family core-pipeline
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
#' print(result)
#' summary(result)
#'
#' @export
der_diagnose <- function(x, ..., tau = 1.2, correct = TRUE) {
  result <- der_compute(x, ...)
  result <- der_classify(result, tau = tau, verbose = FALSE)
  if (correct) {
    result <- der_correct(result)
  }
  result
}
