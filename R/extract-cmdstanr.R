###############################################################################
# extract-cmdstanr.R
# cmdstanr integration for the svyder package
# ---------------------------------------------------------------------------
# Provides extract_draws.CmdStanMCMC() and der_compute.CmdStanMCMC().
# Unlike brms, cmdstanr does NOT store model data in the fit object,
# so all data arguments (y, X, group, weights) must be supplied by the user.
###############################################################################

#' @rdname extract_draws
#'
#' @param pars Character vector of parameter names to extract from the
#'   CmdStanMCMC object. If \code{NULL} (default), extracts all parameters
#'   except \code{lp__} and \code{lp_approx__}.
#'
#' @export
extract_draws.CmdStanMCMC <- function(x, ..., pars = NULL) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Package 'cmdstanr' is required for extract_draws.CmdStanMCMC(). ",
         "Install via: install.packages('cmdstanr', ",
         "repos = c('https://mc-stan.org/r-packages/', getOption('repos')))",
         call. = FALSE)
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for extract_draws.CmdStanMCMC(). ",
         "Install it with: install.packages('posterior')",
         call. = FALSE)
  }

  # Convert to draws matrix using the posterior package
  if (!is.null(pars)) {
    raw_draws <- x$draws(variables = pars)
  } else {
    raw_draws <- x$draws()
  }

  draws_mat <- as.matrix(posterior::as_draws_matrix(raw_draws))

  # Remove diagnostic columns if present
  exclude_cols <- c("lp__", "lp_approx__")
  keep <- !(colnames(draws_mat) %in% exclude_cols)
  draws_mat <- draws_mat[, keep, drop = FALSE]

  list(draws = draws_mat)
}


#' @rdname der_compute
#'
#' @details
#' \strong{CmdStanMCMC method}: CmdStan does not store model data in the fit
#' object, so the user must provide all data arguments: \code{y}, \code{X},
#' \code{group}, \code{weights}, \code{sigma_theta}, and optionally \code{psu},
#' \code{sigma_e} (for gaussian), \code{beta_prior_sd}, and \code{param_types}.
#'
#' The draws are extracted from the CmdStanMCMC object using
#' [extract_draws()] and then passed to the matrix method.
#' The draws matrix must have columns ordered as
#' \code{[beta_1, ..., beta_p, theta_1, ..., theta_J]}.
#'
#' @export
der_compute.CmdStanMCMC <- function(x, ..., y, X, group, weights,
                                     psu = NULL, family = "binomial",
                                     sigma_theta, sigma_e = NULL,
                                     beta_prior_sd = 5,
                                     param_types = NULL,
                                     design = NULL) {

  # Extract draws from the CmdStanMCMC object
  extracted <- extract_draws(x)
  draws_mat <- extracted$draws

  # Delegate to der_compute.matrix with all user-provided arguments
  der_compute.matrix(
    draws_mat, ...,
    y             = y,
    X             = X,
    group         = group,
    weights       = weights,
    psu           = psu,
    family        = family,
    sigma_theta   = sigma_theta,
    sigma_e       = sigma_e,
    beta_prior_sd = beta_prior_sd,
    param_types   = param_types,
    design        = design
  )
}
