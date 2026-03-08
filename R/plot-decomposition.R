###############################################################################
# plot-decomposition.R
# Decomposition plot: DER vs component values
# ---------------------------------------------------------------------------
# Shows the relationship between DER and its constituent factors
# (shrinkage B, design effect DEFF, protection R_k).
###############################################################################


# ggplot2 decomposition plot
.plot_decomposition_gg <- function(x, ...) {

  decomp <- der_decompose(x)

  # Build plot data
  df <- data.frame(
    param      = decomp$param,
    param_type = decomp$param_type,
    der        = decomp$der,
    B_mean     = decomp$B_mean,
    deff_mean  = decomp$deff_mean,
    R_k        = decomp$R_k,
    der_predicted = decomp$der_predicted,
    stringsAsFactors = FALSE
  )

  # Map param_type to display labels
  type_labels <- c(
    "fe_within"  = "Fixed (within)",
    "fe_between" = "Fixed (between)",
    "re"         = "Random effects"
  )
  df$type_label <- type_labels[df$param_type]
  df$type_label[is.na(df$type_label)] <- df$param_type[is.na(df$type_label)]

  # Scatter: x = der_predicted (from decomposition), y = der (observed)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$der_predicted,
                                          y = .data$der,
                                          color = .data$type_label)) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::geom_abline(intercept = 0, slope = 1,
                          linetype = "dashed", color = "gray50") +
    ggplot2::scale_color_manual(
      values = c(
        "Fixed (within)"  = "#E69F00",
        "Fixed (between)" = "#56B4E9",
        "Random effects"  = "#999999"
      ),
      name = "Parameter type"
    ) +
    ggplot2::labs(
      x     = "Predicted DER (from decomposition)",
      y     = "Observed DER",
      title = "DER Decomposition: Observed vs Predicted"
    ) +
    ggplot2::theme_minimal()

  # Add text annotation with mean B and mean DEFF
  B_mean    <- mean(x$B_j)
  deff_mean <- mean(x$deff_j)

  annot_text <- sprintf("Mean B = %.3f\nMean DEFF = %.3f",
                          B_mean, deff_mean)

  p <- p + ggplot2::annotate("text", x = -Inf, y = Inf,
                              label = annot_text,
                              hjust = -0.1, vjust = 1.2,
                              size = 3.5, color = "gray30")

  p
}


# Base R decomposition plot
.plot_decomposition_base <- function(x, ...) {

  decomp <- der_decompose(x)

  df <- data.frame(
    param      = decomp$param,
    param_type = decomp$param_type,
    der        = decomp$der,
    der_predicted = decomp$der_predicted,
    stringsAsFactors = FALSE
  )

  # Colors by param_type
  type_colors <- c(
    "fe_within"  = "#E69F00",
    "fe_between" = "#56B4E9",
    "re"         = "#999999"
  )
  pt_colors <- type_colors[df$param_type]

  xy_range <- range(c(df$der, df$der_predicted))
  xy_lim   <- c(xy_range[1] * 0.9, xy_range[2] * 1.1)

  oldpar <- graphics::par(mar = c(4, 4, 3, 2))
  on.exit(graphics::par(oldpar), add = TRUE)

  graphics::plot(df$der_predicted, df$der,
                  pch = 19, col = pt_colors,
                  xlim = xy_lim, ylim = xy_lim,
                  xlab = "Predicted DER (from decomposition)",
                  ylab = "Observed DER",
                  main = "DER Decomposition: Observed vs Predicted", ...)

  # 1:1 line
  graphics::abline(a = 0, b = 1, lty = 2, col = "gray50")

  # Legend
  types_present <- unique(df$param_type)
  type_labels <- c(
    "fe_within"  = "Fixed (within)",
    "fe_between" = "Fixed (between)",
    "re"         = "Random effects"
  )
  graphics::legend("topleft",
                    legend = type_labels[types_present],
                    col    = type_colors[types_present],
                    pch    = 19, cex = 0.8, bty = "n")

  # Annotation
  B_mean    <- mean(x$B_j)
  deff_mean <- mean(x$deff_j)
  graphics::mtext(sprintf("Mean B = %.3f | Mean DEFF = %.3f",
                            B_mean, deff_mean),
                    side = 3, line = 0, cex = 0.8, col = "gray30")

  invisible(x)
}
