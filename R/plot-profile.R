###############################################################################
# plot-profile.R
# Profile plot: parameters on y-axis, DER on x-axis, tier coloring
# ---------------------------------------------------------------------------
# Both ggplot2 and base R implementations.
###############################################################################

# Colorblind-friendly palette for tiers
# Orange: flagged (I-a CORRECT)
# Blue:   protected between (I-b)
# Gray:   protected random effects (II)
.tier_colors <- c(
  "I-a"  = "#E69F00",
  "I-b"  = "#56B4E9",
  "II"   = "#999999"
)

# Helper: build the profile data.frame used by both ggplot2 and base R
.build_profile_data <- function(x) {
  has_class <- "flagged" %in% names(x$classification)
  d <- length(x$der)

  if (has_class) {
    df <- data.frame(
      param  = x$classification$param_name,
      der    = x$classification$der,
      tier   = x$classification$tier,
      action = x$classification$action,
      stringsAsFactors = FALSE
    )
  } else {
    # No classification: assign all to "II" tier with retain action
    df <- data.frame(
      param  = x$params,
      der    = as.numeric(x$der),
      tier   = rep("II", d),
      action = rep("retain", d),
      stringsAsFactors = FALSE
    )
  }

  # Sort by DER value (ascending)
  df <- df[order(df$der), ]
  df$param <- factor(df$param, levels = df$param)
  df
}


# ggplot2 profile plot
.plot_profile_gg <- function(x, tau = NULL, ...) {

  df <- .build_profile_data(x)

  # Determine tau for threshold line
  if (is.null(tau)) {
    tau <- if (!is.na(x$tau)) x$tau else NULL
  }

  # Determine x-axis scale: use log10 if range spans > 100x
  der_range <- range(df$der)
  use_log <- (der_range[2] / max(der_range[1], 1e-6)) > 100

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$der, y = .data$param,
                                          color = .data$tier)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(
      values = .tier_colors,
      name   = "Tier",
      labels = c("I-a" = "I-a Survey-dominated",
                  "I-b" = "I-b Protected (between)",
                  "II"  = "II Protected (RE)")
    ) +
    ggplot2::labs(
      x     = "Design Effect Ratio (DER)",
      y     = NULL,
      title = "DER Profile Plot"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "gray90"),
      panel.grid.minor   = ggplot2::element_blank()
    )

  # Reference line at DER = 1.0
  p <- p + ggplot2::geom_vline(xintercept = 1.0, linetype = "solid",
                                color = "gray60", linewidth = 0.5)

  # Threshold line

  if (!is.null(tau)) {
    p <- p + ggplot2::geom_vline(xintercept = tau, linetype = "dashed",
                                  color = "red", linewidth = 0.7)
  }

  if (use_log) {
    p <- p + ggplot2::scale_x_log10()
  }

  p
}


# Base R profile plot
.plot_profile_base <- function(x, tau = NULL, ...) {

  df <- .build_profile_data(x)

  if (is.null(tau)) {
    tau <- if (!is.na(x$tau)) x$tau else NULL
  }

  # Colors for each point
  pt_colors <- .tier_colors[df$tier]

  # Determine x-axis range
  xrange <- range(df$der)
  if (!is.null(tau)) {
    xrange <- range(c(xrange, tau))
  }
  xrange <- c(xrange[1] * 0.9, xrange[2] * 1.1)

  n <- nrow(df)
  ypos <- seq_len(n)

  oldpar <- graphics::par(mar = c(4, 10, 3, 2))
  on.exit(graphics::par(oldpar), add = TRUE)

  graphics::plot(df$der, ypos, pch = 19, col = pt_colors,
                  xlim = xrange, ylim = c(0.5, n + 0.5),
                  xlab = "Design Effect Ratio (DER)", ylab = "",
                  yaxt = "n", main = "DER Profile Plot", ...)

  graphics::axis(2, at = ypos, labels = as.character(df$param),
                  las = 1, cex.axis = 0.7)

  # Reference line at DER = 1.0
  graphics::abline(v = 1.0, lty = 1, col = "gray60")

  # Threshold line
  if (!is.null(tau)) {
    graphics::abline(v = tau, lty = 2, col = "red", lwd = 2)
  }

  # Legend
  tiers_present <- unique(df$tier)
  tier_labels <- c("I-a" = "I-a Survey-dominated",
                    "I-b" = "I-b Protected (between)",
                    "II"  = "II Protected (RE)")
  graphics::legend("bottomright",
                    legend = tier_labels[tiers_present],
                    col    = .tier_colors[tiers_present],
                    pch    = 19, cex = 0.8, bty = "n")

  invisible(x)
}
