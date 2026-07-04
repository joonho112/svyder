###############################################################################
# der_decompose.R
# Decompose DER into constituent factors
# ---------------------------------------------------------------------------
# Breaks each parameter's DER into design effect, shrinkage, protection
# factor R_k, and finite-J correction kappa -- using STRUCTURAL quantities
# computed from the design matrix and working weights, never back-solved
# from the observed DER (a back-solved R_k makes der_predicted == der by
# construction and validates nothing).
###############################################################################

#' Decompose DER into Components
#'
#' Decomposes each parameter's DER into its constituent factors. For fixed
#' effects, the protection factor \eqn{R_k} is computed structurally from
#' the working-weighted within/between split of covariate \eqn{k}'s
#' identifying variation:
#' \deqn{R_k = \frac{\sum_j B_j\, a_j \bar x_{jk}^2}
#'                  {\sum_j (W_{jkk} + a_j \bar x_{jk}^2)},}
#' where \eqn{a_j} is the working information of group \eqn{j},
#' \eqn{\bar x_{jk}} the working-weighted group mean of covariate \eqn{k},
#' and \eqn{W_{jkk}} its within-group working sum of squares. A pure
#' within-group covariate has \eqn{R_k = 0}; a pure between-group covariate
#' has \eqn{R_k} equal to the information-weighted average shrinkage factor.
#' The prediction is \eqn{DER_k \approx DEFF \cdot (1 - R_k)} with the
#' global Kish DEFF.
#'
#' For random effects the per-group prediction
#' \eqn{DER_j \approx B_j \cdot DEFF_j \cdot \kappa_j} is used (matching
#' [der_theorem_check()]), not a pooled mean.
#'
#' Because \code{der_predicted} is structural, comparing it with the
#' observed \code{der} is a genuine check of the decomposition theorems'
#' mechanism, within their stated model class.
#'
#' The prediction is row-specific: \code{der_predicted} uses Theorem 1
#' (\eqn{\mathrm{DEFF}\cdot(1 - R_k)}) for fixed-effect rows and Theorem 2
#' (\eqn{B_j\cdot\mathrm{DEFF}_j\cdot\kappa_j}) for random-effect rows.
#' Consequently the component columns are populated per row type:
#' \code{B_used} and \code{kappa} are \code{NA} for fixed effects (whose
#' prediction uses only DEFF and \eqn{R_k}), and \code{R_k} is \code{NA} for
#' random effects (whose prediction uses only \eqn{B_j}, \eqn{\mathrm{DEFF}_j}
#' and \eqn{\kappa_j}).
#'
#' @param x A \code{svyder} object from svyder >= 0.2.0 (the stored
#'   \code{$data} slots are required).
#'
#' @return A \code{data.frame} with one row per parameter and columns:
#'   \describe{
#'     \item{param}{Parameter name.}
#'     \item{param_type}{\code{"fe_within"}, \code{"fe_between"}, or
#'       \code{"re"}.}
#'     \item{der}{Observed Design Effect Ratio.}
#'     \item{deff_used}{Design effect entering the prediction: the global
#'       Kish DEFF for fixed effects, the per-group \eqn{\mathrm{DEFF}_j} for
#'       random effects.}
#'     \item{B_used}{Shrinkage factor \eqn{B_j} (random-effect rows only;
#'       \code{NA} for fixed effects).}
#'     \item{R_k}{Structural between-group share of covariate \eqn{k}'s
#'       identifying variation (fixed-effect rows only; \code{NA} for random
#'       effects).}
#'     \item{kappa}{Finite-group correction \eqn{\kappa_j} (random-effect
#'       rows only; \code{NA} for fixed effects).}
#'     \item{der_predicted}{Theorem-based prediction: \eqn{\mathrm{DEFF}\,
#'       (1 - R_k)} for fixed effects (Theorem 1),
#'       \eqn{B_j\,\mathrm{DEFF}_j\,\kappa_j} for random effects (Theorem 2).}
#'   }
#'
#' @references
#' Lee, J., Williams, M. R., & Savitsky, T. D. (2026). Design Effect Ratios
#' for Bayesian Survey Models: A Diagnostic Framework for Identifying
#' Survey-Sensitive Parameters. \emph{Journal of Survey Statistics and
#' Methodology}. Submitted.
#'
#' Kish, L. (1965). \emph{Survey Sampling}. John Wiley & Sons.
#'
#' Fay, R. E., & Herriot, R. A. (1979). Estimates of income for small places:
#' An application of James--Stein procedures to census data. \emph{Journal of
#' the American Statistical Association}, 74(366), 269--277.
#' \doi{10.1080/01621459.1979.10482505}
#'
#' @seealso [der_theorem_check()] for theorem-level comparison,
#'   [plot.svyder()] with \code{type = "decomposition"} for visualization.
#' @family analysis
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
#' decomp <- der_decompose(result)
#' head(decomp)
#'
#' @export
der_decompose <- function(x) {

  stopifnot(is.svyder(x))

  if (is.null(x$data)) {
    stop("This svyder object has no stored data slots (created by ",
         "svyder < 0.2.0). Re-run der_compute() to use the structural ",
         "decomposition.", call. = FALSE)
  }

  der         <- x$der
  d           <- length(der)
  param_names <- x$params
  param_types <- x$classification$param_type

  J      <- x$n_groups
  deff_j <- x$deff_j
  B_j    <- x$B_j

  X     <- x$data$X
  group <- x$data$group
  v     <- x$data$v
  w     <- x$data$weights
  p     <- ncol(X)

  # --- Global Kish DEFF (weights-only) ---
  deff_global <- length(w) * sum(w^2) / sum(w)^2

  # --- Structural R_k per fixed-effect column ---
  a_j <- as.numeric(tapply(v, group, sum))                    # length J
  R_k <- numeric(p)
  for (k in seq_len(p)) {
    xbar_jk <- as.numeric(tapply(v * X[, k], group, sum)) / a_j
    between_k <- a_j * xbar_jk^2
    total_k   <- sum(v * X[, k]^2)                            # = W_jkk + between
    R_k[k] <- sum(B_j * between_k) / total_k
  }

  # --- Per-group kappa ---
  kappa_j <- .compute_kappa_j(B_j, J)

  # --- Assemble ---
  out_param     <- character(d)
  out_type      <- character(d)
  out_der       <- numeric(d)
  out_deff      <- numeric(d)
  out_B         <- numeric(d)
  out_Rk        <- numeric(d)
  out_kappa     <- numeric(d)
  out_predicted <- numeric(d)

  n_fe <- sum(param_types != "re")

  for (i in seq_len(d)) {
    out_param[i] <- param_names[i]
    out_type[i]  <- param_types[i]
    out_der[i]   <- der[i]

    if (param_types[i] %in% c("fe_within", "fe_between")) {
      out_deff[i]      <- deff_global
      out_B[i]         <- NA_real_
      out_Rk[i]        <- R_k[i]
      out_kappa[i]     <- NA_real_
      out_predicted[i] <- deff_global * (1 - R_k[i])
    } else {
      j_idx <- i - n_fe
      out_deff[i]      <- deff_j[j_idx]
      out_B[i]         <- B_j[j_idx]
      out_Rk[i]        <- NA_real_
      out_kappa[i]     <- kappa_j[j_idx]
      out_predicted[i] <- B_j[j_idx] * deff_j[j_idx] * kappa_j[j_idx]
    }
  }

  data.frame(
    param         = out_param,
    param_type    = out_type,
    der           = out_der,
    deff_used     = out_deff,
    B_used        = out_B,
    R_k           = out_Rk,
    kappa         = out_kappa,
    der_predicted = out_predicted,
    stringsAsFactors = FALSE
  )
}
