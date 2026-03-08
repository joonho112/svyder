###############################################################################
# plot.R
# Main plot dispatcher for svyder objects
# ---------------------------------------------------------------------------
# Dispatches to profile, decomposition, or comparison plots.
# Each type has both ggplot2 and base R implementations.
###############################################################################

#' Plot DER Diagnostic Results
#'
#' Creates diagnostic visualizations for Design Effect Ratio results.
#' Three plot types are available:
#' \describe{
#'   \item{profile}{DER values by parameter with tier-based coloring and
#'     threshold line. Default plot type.}
#'   \item{decomposition}{Scatter of DER against decomposition components
#'     (shrinkage factor B, design effect DEFF).}
#'   \item{comparison}{Side-by-side naive vs corrected credible intervals.
#'     Requires [der_correct()] to have been applied.}
#' }
#'
#' If \pkg{ggplot2} is installed, returns a ggplot object. Otherwise,
#' falls back to base R graphics.
#'
#' @param x A \code{svyder} object.
#' @param type Character; plot type: \code{"profile"} (default),
#'   \code{"decomposition"}, or \code{"comparison"}.
#' @param ... Additional arguments passed to the underlying plot function.
#'
#' @return If \pkg{ggplot2} is available, a \code{ggplot} object (invisibly).
#'   Otherwise, \code{x} is returned invisibly.
#'
#' @seealso [autoplot.svyder()] for explicit ggplot2 usage.
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
#'   plot(result, type = "profile")
#' }
#'
#' @export
plot.svyder <- function(x, type = c("profile", "decomposition", "comparison"), ...) {
  type <- match.arg(type)
  has_gg <- requireNamespace("ggplot2", quietly = TRUE)

  switch(type,
    profile       = if (has_gg) .plot_profile_gg(x, ...) else .plot_profile_base(x, ...),
    decomposition = if (has_gg) .plot_decomposition_gg(x, ...) else .plot_decomposition_base(x, ...),
    comparison    = if (has_gg) .plot_comparison_gg(x, ...) else .plot_comparison_base(x, ...)
  )
}
