###############################################################################
# der_theorem_check.R
# Verify theoretical DER predictions against empirical values
# ---------------------------------------------------------------------------
# Compares empirical DER values (from the sandwich computation) against
# theoretical predictions from the decomposition theorems.
###############################################################################

#' Verify Theoretical DER Predictions
#'
#' Compares empirical DER values against theoretical predictions
#' from the decomposition theorems. This function checks how well
#' the closed-form approximations match the numerically computed DER.
#'
#' For fixed effects (Theorem 1):
#' \itemize{
#'   \item \code{fe_within}: DER \eqn{\approx}{~} DEFF.
#'   \item \code{fe_between}: DER \eqn{\approx}{~} DEFF \eqn{\cdot}{*} (1 - B).
#' }
#'
#' For random effects (Theorem 2):
#' \itemize{
#'   \item DER \eqn{\approx}{~} B \eqn{\cdot}{*} DEFF \eqn{\cdot}{*} kappa(J).
#' }
#'
#' Also checks the conservation law (Corollary 5) when applicable:
#' DER_mu + DER_theta_cond \eqn{\approx}{~} DEFF (balanced intercept-only case).
#'
#' @param x A \code{svyder} object.
#'
#' @return A \code{data.frame} with columns: \code{param}, \code{param_type},
#'   \code{der_empirical}, \code{der_theorem1} (for FE), \code{der_theorem2}
#'   (for RE), \code{relative_error}, \code{theorem_used}. If the conservation
#'   law is applicable, the result has a \code{"conservation_law"} attribute.
#'
#' @seealso [der_decompose()] for the full decomposition.
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
#' thm <- der_theorem_check(result)
#' head(thm)
#'
#' @export
der_theorem_check <- function(x) {


  stopifnot(is.svyder(x))

  der         <- x$der
  d           <- length(der)
  param_names <- x$params
  param_types <- x$classification$param_type

  J      <- x$n_groups
  deff_j <- x$deff_j
  B_j    <- x$B_j

  deff_mean  <- mean(deff_j)
  B_mean     <- mean(B_j)
  kappa_j    <- .compute_kappa_j(B_j, J)
  kappa_mean <- mean(kappa_j)

  # --- Build output vectors ---
  out_param     <- character(d)
  out_type      <- character(d)
  out_empirical <- numeric(d)
  out_thm1      <- numeric(d)
  out_thm2      <- numeric(d)
  out_relerr    <- numeric(d)
  out_theorem   <- character(d)

  for (i in seq_len(d)) {
    out_param[i]     <- param_names[i]
    out_type[i]      <- param_types[i]
    out_empirical[i] <- der[i]

    pt <- param_types[i]

    if (pt == "fe_within") {
      # Theorem 1 (within): DER ~ DEFF, R_k ~ 0
      predicted     <- deff_mean
      out_thm1[i]   <- predicted
      out_thm2[i]   <- NA_real_
      out_theorem[i] <- "Theorem 1 (within)"

    } else if (pt == "fe_between") {
      # Theorem 1 (between): DER ~ DEFF * (1 - B)
      predicted     <- deff_mean * (1 - B_mean)
      out_thm1[i]   <- predicted
      out_thm2[i]   <- NA_real_
      out_theorem[i] <- "Theorem 1 (between)"

    } else {
      # Theorem 2 (RE): DER ~ B * DEFF * kappa(J)
      # Use per-group values for the j-th random effect
      j_idx <- i - sum(param_types[seq_len(i)] != "re") + sum(param_types == "re" & seq_along(param_types) <= i) - sum(param_types[seq_len(i)] == "re") + 1L

      # Simpler approach: random effects are at positions after all FE
      n_fe <- sum(param_types != "re")
      j_idx <- i - n_fe

      if (j_idx >= 1L && j_idx <= J) {
        predicted <- B_j[j_idx] * deff_j[j_idx] * kappa_j[j_idx]
      } else {
        # Fallback to mean values
        predicted <- B_mean * deff_mean * kappa_mean
      }

      out_thm1[i]    <- NA_real_
      out_thm2[i]    <- predicted
      out_theorem[i] <- "Theorem 2 (RE)"
    }

    # Relative error: |empirical - predicted| / |empirical|
    predicted_val <- if (!is.na(out_thm1[i])) out_thm1[i] else out_thm2[i]
    if (abs(out_empirical[i]) > .Machine$double.eps) {
      out_relerr[i] <- abs(out_empirical[i] - predicted_val) / abs(out_empirical[i])
    } else {
      out_relerr[i] <- NA_real_
    }
  }

  result <- data.frame(
    param         = out_param,
    param_type    = out_type,
    der_empirical = out_empirical,
    der_theorem1  = out_thm1,
    der_theorem2  = out_thm2,
    relative_error = out_relerr,
    theorem_used  = out_theorem,
    stringsAsFactors = FALSE
  )

  # --- Conservation law check (Corollary 5) ---
  # DER_mu + DER_theta_cond ~ DEFF

  # Only applicable when there is an intercept (fe_between) and RE
  has_between <- any(param_types == "fe_between")
  has_re      <- any(param_types == "re")

  if (has_between && has_re) {
    # Use first fe_between as the intercept proxy
    idx_between <- which(param_types == "fe_between")[1]
    der_mu      <- der[idx_between]

    # Mean DER of random effects as proxy for conditional DER
    idx_re      <- which(param_types == "re")
    der_theta_mean <- mean(der[idx_re])

    conservation_sum   <- der_mu + der_theta_mean
    conservation_error <- abs(conservation_sum - deff_mean) / deff_mean

    attr(result, "conservation_law") <- list(
      der_mu             = der_mu,
      der_theta_mean     = der_theta_mean,
      conservation_sum   = conservation_sum,
      deff_mean          = deff_mean,
      relative_error     = conservation_error
    )
  }

  result
}
