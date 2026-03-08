###############################################################################
# der_classify.R
# Classify parameters by design sensitivity using DER thresholds
# ---------------------------------------------------------------------------
# Assigns each parameter to a tier and flags those exceeding the threshold.
###############################################################################

#' Classify Parameters by Design Sensitivity
#'
#' Assigns each parameter to a design-sensitivity tier and flags those
#' whose DER exceeds the threshold \code{tau}. The three-tier classification:
#' \itemize{
#'   \item \strong{Tier I-a} (\code{fe_within}): Survey-dominated parameters.
#'   \item \strong{Tier I-b} (\code{fe_between}): Protected between-cluster parameters.
#'   \item \strong{Tier II} (\code{re}): Protected random effects.
#' }
#'
#' Parameters with DER > \code{tau} (strict inequality) are flagged for
#' correction, regardless of tier.
#'
#' @param x A \code{svyder} object from [der_compute()].
#' @param tau Threshold (default 1.2). Parameters with DER > tau are flagged.
#' @param verbose Print classification summary (default \code{TRUE}).
#'
#' @return A \code{svyder} object with updated \code{classification} and
#'   \code{tau} fields.
#'
#' @seealso [der_compute()] for computing DER values, [der_correct()] for
#'   applying corrections, [der_diagnose()] for the all-in-one pipeline.
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
#' result <- der_classify(result, tau = 1.2)
#'
#' @export
der_classify <- function(x, tau = 1.2, verbose = TRUE) {

  if (!is.svyder(x)) {
    stop("'x' must be a svyder object from der_compute().", call. = FALSE)
  }

  stopifnot(is.numeric(tau), length(tau) == 1L, is.finite(tau), tau > 0)

  der         <- x$der
  d           <- length(der)
  param_names <- x$params
  param_types <- x$classification$param_type

  # --- Build tier assignments ---
  tier       <- character(d)
  tier_label <- character(d)
  action     <- character(d)
  flagged    <- logical(d)

  for (i in seq_len(d)) {
    pt <- param_types[i]

    if (pt == "fe_within") {
      tier[i]       <- "I-a"
      tier_label[i] <- "Survey-dominated"
    } else if (pt == "fe_between") {
      tier[i]       <- "I-b"
      tier_label[i] <- "Protected (between)"
    } else {
      tier[i]       <- "II"
      tier_label[i] <- "Protected (random effects)"
    }

    # Strict inequality: DER exactly at tau is NOT flagged
    if (der[i] > tau) {
      action[i]  <- "CORRECT"
      flagged[i] <- TRUE
    } else {
      action[i]  <- "retain"
      flagged[i] <- FALSE
    }
  }

  classification <- data.frame(
    param_name = param_names,
    param_type = param_types,
    der        = as.numeric(der),
    tier       = tier,
    tier_label = tier_label,
    flagged    = flagged,
    action     = action,
    stringsAsFactors = FALSE
  )

  x$classification <- classification
  x$tau            <- tau

  # --- Print summary if verbose ---
  if (verbose) {
    n_flagged <- sum(flagged)
    cat(sprintf("DER Classification (tau = %.2f)\n", tau))
    cat(sprintf("  Total parameters: %d\n", d))
    cat(sprintf("  Flagged: %d (%s)\n", n_flagged, .format_pct(n_flagged / d)))

    if (n_flagged > 0) {
      flagged_df <- classification[flagged, , drop = FALSE]
      cat("  Flagged parameters:\n")
      for (row_i in seq_len(nrow(flagged_df))) {
        cat(sprintf("    %s: DER = %.3f [%s] -> %s\n",
                    flagged_df$param_name[row_i],
                    flagged_df$der[row_i],
                    flagged_df$tier[row_i],
                    flagged_df$action[row_i]))
      }
    }
  }

  x
}
