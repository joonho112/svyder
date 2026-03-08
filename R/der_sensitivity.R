###############################################################################
# der_sensitivity.R
# Sensitivity analysis across threshold values
# ---------------------------------------------------------------------------
# Evaluates how the number of flagged parameters changes across a range
# of threshold values tau.
###############################################################################

#' Sensitivity Analysis Across Threshold Values
#'
#' Evaluates how the number of flagged parameters changes across
#' a range of threshold values \code{tau}. Useful for assessing the
#' robustness of classification results to the choice of threshold.
#'
#' @param x A \code{svyder} object.
#' @param tau_range Numeric vector of threshold values to evaluate.
#'   Default: \code{seq(0.8, 2.0, by = 0.1)}.
#'
#' @return A \code{data.frame} with columns: \code{tau}, \code{n_flagged},
#'   \code{pct_flagged}, and \code{flagged_params} (a list-column of
#'   character vectors naming the flagged parameters at each threshold).
#'
#' @seealso [der_classify()] for classification at a single threshold.
#' @family analysis
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
#' sens <- der_sensitivity(result)
#' sens[, c("tau", "n_flagged", "pct_flagged")]
#'
#' @export
der_sensitivity <- function(x, tau_range = seq(0.8, 2.0, by = 0.1)) {

  stopifnot(is.svyder(x))
  stopifnot(is.numeric(tau_range), length(tau_range) >= 1L)

  der         <- x$der
  d           <- length(der)
  param_names <- x$params
  n_tau       <- length(tau_range)

  # Sort tau_range to ensure monotone output
  tau_sorted <- sort(tau_range)

  out_tau     <- numeric(n_tau)
  out_nflag   <- integer(n_tau)
  out_pct     <- numeric(n_tau)
  out_params  <- vector("list", n_tau)

  for (k in seq_len(n_tau)) {
    tau_k <- tau_sorted[k]
    # Strict inequality: DER > tau flags the parameter
    flagged_k <- der > tau_k
    flagged_names <- param_names[flagged_k]

    out_tau[k]    <- tau_k
    out_nflag[k]  <- sum(flagged_k)
    out_pct[k]    <- sum(flagged_k) / d
    out_params[[k]] <- flagged_names
  }

  data.frame(
    tau            = out_tau,
    n_flagged      = out_nflag,
    pct_flagged    = out_pct,
    flagged_params = I(out_params),
    stringsAsFactors = FALSE
  )
}
