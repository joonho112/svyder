###############################################################################
# plot-comparison.R
# CI comparison plot: naive vs corrected intervals
# ---------------------------------------------------------------------------
# Shows side-by-side credible intervals before and after DER correction.
# Only meaningful after der_correct() has been applied.
###############################################################################


# ggplot2 comparison plot
.plot_comparison_gg <- function(x, params = NULL, prob = 0.95, ...) {

  # Validate that correction has been applied
  if (is.null(x$corrected_draws)) {
    stop("Comparison plot requires corrected draws. Run der_correct() first.",
         call. = FALSE)
  }

  d <- length(x$der)

  # Determine which parameters to show
  if (is.null(params)) {
    # Show only flagged parameters by default
    has_class <- "flagged" %in% names(x$classification)
    if (has_class) {
      show_idx <- which(x$classification$flagged)
    } else {
      show_idx <- seq_len(d)
    }
    # Fallback: if none flagged, show all
    if (length(show_idx) == 0) {
      show_idx <- seq_len(d)
    }
  } else {
    # Match user-specified parameter names
    show_idx <- match(params, x$params)
    show_idx <- show_idx[!is.na(show_idx)]
    if (length(show_idx) == 0) {
      stop("No matching parameters found.", call. = FALSE)
    }
  }

  # Compute CIs for original and corrected draws
  ci_naive     <- .compute_ci(x$original_draws[, show_idx, drop = FALSE], prob = prob)
  ci_corrected <- .compute_ci(x$corrected_draws[, show_idx, drop = FALSE], prob = prob)

  mean_naive     <- colMeans(x$original_draws[, show_idx, drop = FALSE])
  mean_corrected <- colMeans(x$corrected_draws[, show_idx, drop = FALSE])

  param_labels <- x$params[show_idx]
  n_show <- length(show_idx)

  # Build data.frame for ggplot
  df <- data.frame(
    param    = rep(param_labels, 2),
    estimate = c(mean_naive, mean_corrected),
    lower    = c(ci_naive[, "lower"], ci_corrected[, "lower"]),
    upper    = c(ci_naive[, "upper"], ci_corrected[, "upper"]),
    type     = rep(c("Naive", "Corrected"), each = n_show),
    stringsAsFactors = FALSE
  )

  # Order parameters by DER value
  der_order <- x$der[show_idx]
  param_ordered <- param_labels[order(der_order)]
  df$param <- factor(df$param, levels = param_ordered)

  # Slight vertical dodge so intervals don't overlap
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$estimate, y = .data$param,
                                          xmin = .data$lower, xmax = .data$upper,
                                          color = .data$type)) +
    ggplot2::geom_pointrange(
      position = ggplot2::position_dodge(width = 0.4),
      size = 0.5
    ) +
    ggplot2::scale_color_manual(
      values = c("Naive" = "#999999", "Corrected" = "#E69F00"),
      name   = "Interval"
    ) +
    ggplot2::labs(
      x     = "Estimate",
      y     = NULL,
      title = sprintf("Credible Interval Comparison (%s%%)", prob * 100)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "gray90"),
      panel.grid.minor   = ggplot2::element_blank()
    )

  p
}


# Base R comparison plot
.plot_comparison_base <- function(x, params = NULL, prob = 0.95, ...) {

  if (is.null(x$corrected_draws)) {
    stop("Comparison plot requires corrected draws. Run der_correct() first.",
         call. = FALSE)
  }

  d <- length(x$der)

  # Determine which parameters to show
  if (is.null(params)) {
    has_class <- "flagged" %in% names(x$classification)
    if (has_class) {
      show_idx <- which(x$classification$flagged)
    } else {
      show_idx <- seq_len(d)
    }
    if (length(show_idx) == 0) {
      show_idx <- seq_len(d)
    }
  } else {
    show_idx <- match(params, x$params)
    show_idx <- show_idx[!is.na(show_idx)]
    if (length(show_idx) == 0) {
      stop("No matching parameters found.", call. = FALSE)
    }
  }

  # Compute CIs
  ci_naive     <- .compute_ci(x$original_draws[, show_idx, drop = FALSE], prob = prob)
  ci_corrected <- .compute_ci(x$corrected_draws[, show_idx, drop = FALSE], prob = prob)

  mean_naive     <- colMeans(x$original_draws[, show_idx, drop = FALSE])
  mean_corrected <- colMeans(x$corrected_draws[, show_idx, drop = FALSE])

  param_labels <- x$params[show_idx]
  n_show <- length(show_idx)

  # Sort by DER
  der_order <- order(x$der[show_idx])
  param_labels   <- param_labels[der_order]
  mean_naive     <- mean_naive[der_order]
  mean_corrected <- mean_corrected[der_order]
  ci_naive       <- ci_naive[der_order, , drop = FALSE]
  ci_corrected   <- ci_corrected[der_order, , drop = FALSE]

  # X-axis range
  xrange <- range(c(ci_naive, ci_corrected))
  xrange <- c(xrange[1] - 0.1 * diff(xrange), xrange[2] + 0.1 * diff(xrange))

  ypos <- seq_len(n_show)
  dodge <- 0.15

  oldpar <- graphics::par(mar = c(4, 10, 3, 2))
  on.exit(graphics::par(oldpar), add = TRUE)

  graphics::plot(NULL, xlim = xrange, ylim = c(0.5, n_show + 0.5),
                  xlab = "Estimate", ylab = "",
                  yaxt = "n",
                  main = sprintf("Credible Interval Comparison (%s%%)", prob * 100),
                  ...)

  graphics::axis(2, at = ypos, labels = param_labels, las = 1, cex.axis = 0.7)

  # Naive intervals (gray)
  for (i in seq_len(n_show)) {
    y_i <- ypos[i] - dodge
    graphics::segments(ci_naive[i, "lower"], y_i, ci_naive[i, "upper"], y_i,
                        col = "#999999", lwd = 2)
    graphics::points(mean_naive[i], y_i, pch = 19, col = "#999999", cex = 1.2)
  }

  # Corrected intervals (orange)
  for (i in seq_len(n_show)) {
    y_i <- ypos[i] + dodge
    graphics::segments(ci_corrected[i, "lower"], y_i, ci_corrected[i, "upper"], y_i,
                        col = "#E69F00", lwd = 2)
    graphics::points(mean_corrected[i], y_i, pch = 19, col = "#E69F00", cex = 1.2)
  }

  graphics::legend("bottomright",
                    legend = c("Naive", "Corrected"),
                    col    = c("#999999", "#E69F00"),
                    lwd    = 2, pch = 19, cex = 0.8, bty = "n")

  invisible(x)
}
