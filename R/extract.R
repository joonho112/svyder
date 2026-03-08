###############################################################################
# extract.R
# S3 generics and basic methods for extracting posterior draws and survey
# design information from model objects
# ---------------------------------------------------------------------------
# Generics: extract_draws(), extract_design()
# Basic methods: matrix, default
# Backend-specific methods live in extract-brms.R, extract-cmdstanr.R, etc.
###############################################################################

#' Extract Posterior Draws from Model Objects
#'
#' S3 generic for extracting posterior draws from fitted model objects.
#' Each method returns a standardized list with a draws matrix and optional
#' parameter metadata, which can then be passed to [der_compute()].
#'
#' Methods are provided for \code{matrix} (identity), \code{brmsfit},
#' \code{CmdStanMCMC}, \code{stanreg}, \code{draws_matrix}, \code{draws_df},
#' \code{draws_array}, \code{draws_list}, and \code{draws_rvars}.
#'
#' @param x A fitted model object (matrix, brmsfit, CmdStanMCMC, stanreg,
#'   draws_matrix, or draws_df).
#' @param ... Additional arguments passed to methods.
#'
#' @return A list with components:
#'   \describe{
#'     \item{draws}{Numeric matrix of posterior draws (S x d), where S is the
#'       number of posterior samples and d is the number of parameters.}
#'     \item{param_info}{(Optional) Data frame with parameter metadata
#'       including original names and mapped names.}
#'   }
#'
#' @seealso [extract_design()] for extracting survey design information,
#'   [der_compute()] for computing DER from draws.
#' @family extraction
#'
#' @examples
#' # Matrix method is the identity
#' m <- matrix(rnorm(100), nrow = 20)
#' result <- extract_draws(m)
#' identical(result$draws, m)
#'
#' @export
extract_draws <- function(x, ...) UseMethod("extract_draws")

#' @rdname extract_draws
#' @export
extract_draws.matrix <- function(x, ...) {
  # Identity: the matrix IS the draws
  list(draws = x)
}

#' @rdname extract_draws
#' @export
extract_draws.default <- function(x, ...) {
  stop("extract_draws() does not know how to handle class '",
       class(x)[1], "'.\n",
       "Supported: matrix, brmsfit, CmdStanMCMC, stanreg, draws_matrix, draws_df.",
       call. = FALSE)
}


#' Extract Survey Design Information
#'
#' Extracts weights, cluster identifiers, and strata from survey design
#' objects created by the \pkg{survey} package. The extracted information
#' can be passed directly to [der_compute()].
#'
#' @param design A survey design object (class \code{survey.design2}).
#' @param ... Additional arguments passed to methods.
#'
#' @return A list with components:
#'   \describe{
#'     \item{weights}{Numeric vector of survey weights.}
#'     \item{cluster}{Factor or integer vector of PSU/cluster identifiers.}
#'     \item{strata}{Data frame of strata variables (may be \code{NULL}
#'       for unstratified designs).}
#'   }
#'
#' @seealso [extract_draws()] for extracting posterior draws,
#'   [der_compute()] for computing DER.
#' @family extraction
#'
#' @export
extract_design <- function(design, ...) UseMethod("extract_design")

#' @rdname extract_design
#' @export
extract_design.survey.design2 <- function(design, ...) {
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("Package 'survey' is required for extract_design.survey.design2(). ",
         "Install it with: install.packages('survey')",
         call. = FALSE)
  }
  list(
    weights = stats::weights(design),
    cluster = design$cluster[[1]],
    strata  = design$strata
  )
}

#' @rdname extract_design
#' @export
extract_design.default <- function(design, ...) {
  stop("extract_design() does not support class '", class(design)[1], "'.\n",
       "Supported: survey.design2 (from the survey package).",
       call. = FALSE)
}
