# classes.R
# S3 class infrastructure for svyder objects
# ---------------------------------------------------------------------------
# The svyder class wraps the full output of the DER computation pipeline
# into a structured S3 object with validation.


#' Low-Level Constructor for svyder Objects
#'
#' Creates a svyder object from pre-computed DER pipeline results.
#' This is the internal constructor; users should call [der_compute()]
#' or [der_diagnose()] instead.
#'
#' @param der Named numeric vector of Design Effect Ratios.
#' @param params Character vector of parameter names.
#' @param H_obs Numeric matrix, observed information matrix (d x d).
#' @param J_cluster Numeric matrix, clustered score outer product (d x d).
#' @param V_sand Numeric matrix, sandwich variance (d x d).
#' @param sigma_mcmc Numeric matrix, MCMC posterior covariance (d x d).
#' @param deff_j Numeric vector of cluster-level design effects (length J).
#' @param B_j Numeric vector of shrinkage factors (length J).
#' @param classification Data frame with per-parameter classification results.
#' @param tau Numeric scalar, flagging threshold for DER.
#' @param corrected_draws Numeric matrix, corrected posterior draws (M x d).
#' @param scale_factors Numeric vector of Cholesky scale factors (length d).
#' @param original_draws Numeric matrix, original posterior draws (M x d).
#' @param call The matched call that produced this object.
#' @param family Character string, model family (e.g., "binomial", "gaussian").
#' @param n_obs Integer, number of observations.
#' @param n_groups Integer, number of groups/clusters.
#' @param compute_time Numeric scalar, computation time in seconds.
#'
#' @return A list of class \code{"svyder"}.
#'
#' @keywords internal
new_svyder <- function(der,
                       params,
                       H_obs,
                       J_cluster,
                       V_sand,
                       sigma_mcmc,
                       deff_j,
                       B_j,
                       classification,
                       tau,
                       corrected_draws,
                       scale_factors,
                       original_draws,
                       call,
                       family,
                       n_obs,
                       n_groups,
                       compute_time) {

  structure(
    list(
      der              = der,
      params           = params,
      H_obs            = H_obs,
      J_cluster        = J_cluster,
      V_sand           = V_sand,
      sigma_mcmc       = sigma_mcmc,
      deff_j           = deff_j,
      B_j              = B_j,
      classification   = classification,
      tau              = tau,
      corrected_draws  = corrected_draws,
      scale_factors    = scale_factors,
      original_draws   = original_draws,
      call             = call,
      family           = family,
      n_obs            = n_obs,
      n_groups         = n_groups,
      compute_time     = compute_time
    ),
    class = "svyder"
  )
}


#' Validate a svyder Object
#'
#' Checks structural integrity of a svyder object: matching dimensions,
#' correct types, and internal consistency.
#'
#' @param x A svyder object to validate.
#'
#' @return The input object (invisibly) if valid; otherwise throws an error.
#'
#' @keywords internal
validate_svyder <- function(x) {

  if (!is.svyder(x)) {
    stop("Object is not of class 'svyder'.", call. = FALSE)
  }

  # --- Required fields ---
  required <- c("der", "params", "H_obs", "J_cluster", "V_sand",
                 "sigma_mcmc", "deff_j", "B_j", "classification",
                 "tau", "corrected_draws", "scale_factors",
                 "original_draws", "family", "n_obs", "n_groups")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop("Missing required fields: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  d <- length(x$der)

  # --- DER / params length consistency ---
  if (length(x$params) != d) {
    stop("Length of 'params' (", length(x$params),
         ") does not match length of 'der' (", d, ").",
         call. = FALSE)
  }

  # --- Matrix dimension checks ---
  if (!is.matrix(x$H_obs) || nrow(x$H_obs) != d || ncol(x$H_obs) != d) {
    stop("'H_obs' must be a ", d, " x ", d, " matrix.", call. = FALSE)
  }
  if (!is.matrix(x$J_cluster) || nrow(x$J_cluster) != d || ncol(x$J_cluster) != d) {
    stop("'J_cluster' must be a ", d, " x ", d, " matrix.", call. = FALSE)
  }
  if (!is.matrix(x$V_sand) || nrow(x$V_sand) != d || ncol(x$V_sand) != d) {
    stop("'V_sand' must be a ", d, " x ", d, " matrix.", call. = FALSE)
  }
  if (!is.matrix(x$sigma_mcmc) || nrow(x$sigma_mcmc) != d || ncol(x$sigma_mcmc) != d) {
    stop("'sigma_mcmc' must be a ", d, " x ", d, " matrix.", call. = FALSE)
  }

  # --- Classification data frame ---
  if (!is.data.frame(x$classification)) {
    stop("'classification' must be a data.frame.", call. = FALSE)
  }
  if (nrow(x$classification) != d) {
    stop("Number of rows in 'classification' (", nrow(x$classification),
         ") does not match length of 'der' (", d, ").",
         call. = FALSE)
  }

  # --- Draws dimensions ---
  if (!is.null(x$corrected_draws)) {
    if (!is.matrix(x$corrected_draws) || ncol(x$corrected_draws) != d) {
      stop("'corrected_draws' must be a matrix with ", d, " columns.",
           call. = FALSE)
    }
  }
  if (!is.null(x$original_draws)) {
    if (!is.matrix(x$original_draws) || ncol(x$original_draws) != d) {
      stop("'original_draws' must be a matrix with ", d, " columns.",
           call. = FALSE)
    }
  }

  # --- Scale factors ---
  if (length(x$scale_factors) != d) {
    stop("Length of 'scale_factors' (", length(x$scale_factors),
         ") does not match length of 'der' (", d, ").",
         call. = FALSE)
  }

  # --- Scalar checks ---
  if (!is.numeric(x$tau) || length(x$tau) != 1L || !is.finite(x$tau)) {
    stop("'tau' must be a finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(x$n_obs) || length(x$n_obs) != 1L) {
    stop("'n_obs' must be a numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(x$n_groups) || length(x$n_groups) != 1L) {
    stop("'n_groups' must be a numeric scalar.", call. = FALSE)
  }

  # --- deff_j / B_j consistency with n_groups ---
  J <- x$n_groups
  if (length(x$deff_j) != J) {
    stop("Length of 'deff_j' (", length(x$deff_j),
         ") does not match 'n_groups' (", J, ").",
         call. = FALSE)
  }
  if (length(x$B_j) != J) {
    stop("Length of 'B_j' (", length(x$B_j),
         ") does not match 'n_groups' (", J, ").",
         call. = FALSE)
  }

  invisible(x)
}


#' Test if an Object is a svyder Object
#'
#' Checks whether an R object inherits from class \code{"svyder"}.
#'
#' @param x An R object.
#'
#' @return Logical; \code{TRUE} if \code{x} inherits from class
#'   \code{"svyder"}, \code{FALSE} otherwise.
#'
#' @family svyder-methods
#'
#' @examples
#' is.svyder(list())
#'
#' data(nsece_demo)
#' result <- der_compute(
#'   nsece_demo$draws,
#'   y = nsece_demo$y, X = nsece_demo$X,
#'   group = nsece_demo$group, weights = nsece_demo$weights,
#'   psu = nsece_demo$psu, family = "binomial",
#'   sigma_theta = nsece_demo$sigma_theta,
#'   param_types = nsece_demo$param_types
#' )
#' is.svyder(result)
#'
#' @export
is.svyder <- function(x) {
  inherits(x, "svyder")
}
