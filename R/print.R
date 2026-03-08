###############################################################################
# print.R
# Print and summary methods for svyder objects
###############################################################################

#' Print a svyder Object
#'
#' Displays a concise summary of DER diagnostic results, including the
#' DER range, classification threshold, number of flagged parameters,
#' and correction status.
#'
#' @param x A \code{svyder} object.
#' @param n Maximum number of flagged parameters to display (default 10).
#' @param digits Number of decimal places for DER values (default 3).
#' @param ... Ignored.
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso [summary.svyder()] for detailed classification output,
#'   [tidy.svyder()] for a tidy data frame.
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
#' print(result)
#'
#' @export
print.svyder <- function(x, n = 10, digits = 3, ...) {

  d <- length(x$der)

  cat(sprintf("svyder diagnostic (%d parameters)\n", d))
  cat(sprintf("  Family: %s | N = %d | J = %d\n",
              x$family, x$n_obs, x$n_groups))
  cat(sprintf("  DER range: [%.*f, %.*f]\n",
              digits, min(x$der),
              digits, max(x$der)))

  # --- Check if classification has been done ---
  has_classification <- "flagged" %in% names(x$classification)

  if (has_classification) {
    flagged   <- x$classification$flagged
    n_flagged <- sum(flagged)

    cat(sprintf("  Threshold (tau): %.2f\n", x$tau))
    cat(sprintf("  Flagged: %d / %d (%s)\n",
                n_flagged, d, .format_pct(n_flagged / d)))

    if (n_flagged > 0) {
      cat("\n  Flagged parameters:\n")
      flagged_df <- x$classification[flagged, , drop = FALSE]

      # Limit display
      show_n <- min(n, nrow(flagged_df))
      for (i in seq_len(show_n)) {
        cat(sprintf("    %-20s DER = %.*f  [%s] -> %s\n",
                    flagged_df$param_name[i],
                    digits,
                    flagged_df$der[i],
                    flagged_df$tier[i],
                    flagged_df$action[i]))
      }
      if (nrow(flagged_df) > show_n) {
        cat(sprintf("    ... and %d more\n", nrow(flagged_df) - show_n))
      }
    }

    # --- Correction status ---
    if (!is.null(x$corrected_draws)) {
      n_corrected <- sum(x$scale_factors != 1.0)
      cat(sprintf("\n  Correction applied: %d parameter(s) rescaled\n",
                  n_corrected))
    }
  } else {
    cat("  (not yet classified -- run der_classify())\n")
  }

  if (!is.null(x$compute_time) && is.finite(x$compute_time)) {
    cat(sprintf("  Compute time: %.3f sec\n", x$compute_time))
  }

  invisible(x)
}


#' Summarize a svyder Object
#'
#' Returns a data frame with per-parameter classification details,
#' including tier assignment, DER value, and flagging status.
#'
#' @param object A \code{svyder} object.
#' @param ... Ignored.
#'
#' @return A \code{data.frame} with classification details (printed and
#'   returned invisibly).
#'
#' @seealso [print.svyder()] for a concise summary,
#'   [tidy.svyder()] for a tidy data frame with posterior summaries.
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
#' summary(result)
#'
#' @export
summary.svyder <- function(object, ...) {

  has_classification <- "flagged" %in% names(object$classification)

  if (has_classification) {
    out <- object$classification
  } else {
    # Minimal summary: params + DER values
    out <- data.frame(
      param_name = object$params,
      param_type = object$classification$param_type,
      der        = as.numeric(object$der),
      stringsAsFactors = FALSE
    )
  }

  print(out, row.names = FALSE)
  invisible(out)
}
