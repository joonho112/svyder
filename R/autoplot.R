###############################################################################
# autoplot.R
# ggplot2 autoplot method for svyder objects
# ---------------------------------------------------------------------------
# Provides a tidyverse-native entry point for DER visualization.
###############################################################################

#' Create a ggplot2 Visualization of DER Results
#'
#' Provides the same functionality as [plot.svyder()] but
#' requires \pkg{ggplot2} and always returns a \code{ggplot} object.
#' This is the preferred method when using ggplot2 directly.
#'
#' @param object A \code{svyder} object.
#' @param type Character; plot type: \code{"profile"} (default),
#'   \code{"decomposition"}, or \code{"comparison"}.
#' @param ... Additional arguments passed to the underlying plot function.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso [plot.svyder()] for the generic plot method.
#' @family visualization
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
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ggplot2::autoplot(result, type = "profile")
#' }
#'
#' @export
autoplot.svyder <- function(object, type = c("profile", "decomposition", "comparison"), ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for autoplot(). ",
         "Install it with: install.packages('ggplot2')",
         call. = FALSE)
  }
  type <- match.arg(type)
  switch(type,
    profile       = .plot_profile_gg(object, ...),
    decomposition = .plot_decomposition_gg(object, ...),
    comparison    = .plot_comparison_gg(object, ...)
  )
}
