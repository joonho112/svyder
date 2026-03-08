###############################################################################
# der_decompose.R
# Decompose DER into constituent factors
# ---------------------------------------------------------------------------
# Breaks each parameter's DER into Kish DEFF, shrinkage factor B,
# protection factor R_k, and finite-J correction kappa.
###############################################################################

#' Decompose DER into Components
#'
#' Decomposes each parameter's DER into its constituent factors:
#' Kish DEFF, shrinkage factor B, protection factor R_k, and
#' finite-J correction kappa. This decomposition reveals why each
#' parameter has its observed DER value.
#'
#' For fixed effects:
#' \itemize{
#'   \item \code{fe_within}: DER \eqn{\approx}{~} DEFF \eqn{\cdot}{*} (1 - R_k),
#'     where R_k \eqn{\approx}{~} 0.
#'   \item \code{fe_between}: DER \eqn{\approx}{~} DEFF \eqn{\cdot}{*} (1 - R_k),
#'     where R_k \eqn{\approx}{~} B.
#' }
#'
#' For random effects:
#' \itemize{
#'   \item DER \eqn{\approx}{~} B \eqn{\cdot}{*} DEFF \eqn{\cdot}{*} kappa(J).
#' }
#'
#' @param x A \code{svyder} object.
#'
#' @return A \code{data.frame} with columns: \code{param}, \code{param_type},
#'   \code{der}, \code{deff_mean}, \code{B_mean}, \code{R_k}, \code{kappa},
#'   \code{der_predicted}.
#'
#' @seealso [der_theorem_check()] for verifying theoretical predictions,
#'   [plot.svyder()] with \code{type = "decomposition"} for visualization.
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
#' decomp <- der_decompose(result)
#' head(decomp)
#'
#' @export
der_decompose <- function(x) {

  stopifnot(is.svyder(x))

  der         <- x$der
  d           <- length(der)
  param_names <- x$params
  param_types <- x$classification$param_type

  J      <- x$n_groups
  deff_j <- x$deff_j
  B_j    <- x$B_j

  # Representative summary values

  deff_mean <- mean(deff_j)
  B_mean    <- mean(B_j)
  kappa_j   <- .compute_kappa_j(B_j, J)
  kappa_mean <- mean(kappa_j)

  # Build output row by row
  out_param     <- character(d)
  out_type      <- character(d)
  out_der       <- numeric(d)
  out_deff      <- numeric(d)
  out_B         <- numeric(d)
  out_Rk        <- numeric(d)
  out_kappa     <- numeric(d)
  out_predicted <- numeric(d)

  for (i in seq_len(d)) {
    out_param[i] <- param_names[i]
    out_type[i]  <- param_types[i]
    out_der[i]   <- der[i]
    out_deff[i]  <- deff_mean
    out_B[i]     <- B_mean

    pt <- param_types[i]

    if (pt == "fe_within") {
      # Within-cluster covariate: R_k ~ 0, DER ~ DEFF
      # Back-solve R_k from empirical DER: R_k = 1 - DER / DEFF
      rk <- 1 - der[i] / deff_mean
      # Clamp to [0, 1]
      rk <- max(0, min(1, rk))
      out_Rk[i]        <- rk
      out_kappa[i]     <- NA_real_
      out_predicted[i] <- deff_mean * (1 - rk)

    } else if (pt == "fe_between") {
      # Between-cluster covariate: R_k ~ B, DER ~ DEFF * (1 - B)
      # Back-solve R_k from empirical DER: R_k = 1 - DER / DEFF
      rk <- 1 - der[i] / deff_mean
      rk <- max(0, min(1, rk))
      out_Rk[i]        <- rk
      out_kappa[i]     <- NA_real_
      out_predicted[i] <- deff_mean * (1 - rk)

    } else {
      # Random effects: DER ~ B * DEFF * kappa(J)
      out_Rk[i]        <- NA_real_
      out_kappa[i]     <- kappa_mean
      out_predicted[i] <- B_mean * deff_mean * kappa_mean
    }
  }

  data.frame(
    param         = out_param,
    param_type    = out_type,
    der           = out_der,
    deff_mean     = out_deff,
    B_mean        = out_B,
    R_k           = out_Rk,
    kappa         = out_kappa,
    der_predicted = out_predicted,
    stringsAsFactors = FALSE
  )
}
