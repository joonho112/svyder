###############################################################################
# extract-rstanarm.R
# rstanarm integration for the svyder package
# ---------------------------------------------------------------------------
# Provides extract_draws.stanreg() and der_compute.stanreg().
# rstanarm stores the data, formula, and family in the fit object,
# enabling auto-detection of y, X, group, and family.
# Survey weights and PSU must still be provided by the user.
###############################################################################

#' @rdname extract_draws
#'
#' @param pars Character vector of parameter names to extract from the
#'   stanreg object. If \code{NULL} (default), extracts all fixed effect
#'   and group-level parameters, excluding variance components and
#'   log-posterior.
#'
#' @export
extract_draws.stanreg <- function(x, ..., pars = NULL) {
  if (!requireNamespace("rstanarm", quietly = TRUE)) {
    stop("Package 'rstanarm' is required for extract_draws.stanreg(). ",
         "Install it with: install.packages('rstanarm')",
         call. = FALSE)
  }

  # Extract all posterior draws as a matrix
  draws_mat <- as.matrix(x)

  if (!is.null(pars)) {
    # User-specified parameter selection
    missing_pars <- setdiff(pars, colnames(draws_mat))
    if (length(missing_pars) > 0L) {
      stop("Parameters not found in model: ",
           paste(missing_pars, collapse = ", "),
           call. = FALSE)
    }
    draws_mat <- draws_mat[, pars, drop = FALSE]
    param_info <- data.frame(
      original_name = pars,
      svyder_name   = pars,
      type          = .classify_rstanarm_param(pars),
      stringsAsFactors = FALSE
    )
  } else {
    # Auto-detect parameters
    all_names <- colnames(draws_mat)
    mapped <- .map_rstanarm_params(all_names)

    if (nrow(mapped) == 0L) {
      stop("No fixed or group-level parameters found in stanreg object.",
           call. = FALSE)
    }

    draws_mat <- draws_mat[, mapped$original_name, drop = FALSE]
    colnames(draws_mat) <- mapped$svyder_name
    param_info <- mapped
  }

  list(
    draws      = draws_mat,
    param_info = param_info
  )
}


#' Classify rstanarm parameter names
#'
#' @param names Character vector of rstanarm parameter names.
#' @return Character vector of types: "fixed", "random", or "other".
#' @keywords internal
.classify_rstanarm_param <- function(names) {
  types <- rep("other", length(names))
  # Fixed effects: (Intercept), x1, x2, etc. (no prefix like b_ or r_)
  # Random effects: b[(Intercept) group:level]
  types[grepl("^b\\[", names)] <- "random"
  # Variance components
  types[grepl("^Sigma\\[", names)] <- "variance"
  types[names == "sigma"] <- "variance"
  types[names == "log-posterior"] <- "diagnostic"
  # Everything else that is not variance or diagnostic is fixed
  types[types == "other"] <- "fixed"
  types
}


#' Map rstanarm parameter names to svyder convention
#'
#' Converts rstanarm-style parameter names to the svyder naming convention:
#' \itemize{
#'   \item \code{(Intercept)} -> \code{beta[1]}
#'   \item \code{x1} -> \code{beta[2]}
#'   \item \code{b[(Intercept) group:level1]} -> \code{theta[1]}
#'   \item \code{b[(Intercept) group:level2]} -> \code{theta[2]}
#' }
#'
#' Excludes variance components (\code{Sigma[...]}, \code{sigma}) and
#' log-posterior.
#'
#' @param param_names Character vector of rstanarm parameter names.
#'
#' @return Data frame with columns: original_name, svyder_name, type.
#'
#' @keywords internal
.map_rstanarm_params <- function(param_names) {
  # --- Identify fixed effects ---
  # rstanarm fixed effects are plain names: (Intercept), x1, x2, ...
  # Exclude: b[...] (random), Sigma[...] (variance), sigma, log-posterior
  exclude_pattern <- "^(b\\[|Sigma\\[|sigma$|log-posterior$|mean_PPD$|accept_stat__|"
  exclude_pattern <- paste0(exclude_pattern, "stepsize__|treedepth__|")
  exclude_pattern <- paste0(exclude_pattern, "n_leapfrog__|divergent__|energy__)")
  fe_names <- param_names[!grepl(exclude_pattern, param_names)]

  # --- Identify random effects ---
  # rstanarm random effects: b[(Intercept) group:level]
  re_names <- param_names[grepl("^b\\[", param_names)]

  # --- Build mapping ---
  result <- data.frame(
    original_name = character(0),
    svyder_name   = character(0),
    type          = character(0),
    stringsAsFactors = FALSE
  )

  # Map fixed effects: (Intercept) -> beta[1], x1 -> beta[2], ...
  if (length(fe_names) > 0L) {
    fe_mapped <- paste0("beta[", seq_along(fe_names), "]")
    result <- rbind(result, data.frame(
      original_name = fe_names,
      svyder_name   = fe_mapped,
      type          = "fixed",
      stringsAsFactors = FALSE
    ))
  }

  # Map random effects: b[(Intercept) group:level1] -> theta[1], ...
  if (length(re_names) > 0L) {
    re_mapped <- paste0("theta[", seq_along(re_names), "]")
    result <- rbind(result, data.frame(
      original_name = re_names,
      svyder_name   = re_mapped,
      type          = "random",
      stringsAsFactors = FALSE
    ))
  }

  result
}


#' @rdname der_compute
#'
#' @details
#' \strong{stanreg method}: The rstanarm method auto-detects the response
#' vector \code{y}, design matrix \code{X}, grouping variable \code{group},
#' model \code{family}, random effect SD \code{sigma_theta}, and residual SD
#' \code{sigma_e} (for gaussian) from the fitted model object. The user must
#' provide \code{weights} (survey weights) and optionally \code{psu} (primary
#' sampling unit indicators).
#'
#' @export
der_compute.stanreg <- function(x, ..., weights, psu = NULL,
                                 sigma_theta = NULL,
                                 sigma_e = NULL,
                                 beta_prior_sd = 5,
                                 param_types = NULL,
                                 design = NULL) {
  if (!requireNamespace("rstanarm", quietly = TRUE)) {
    stop("Package 'rstanarm' is required for der_compute.stanreg(). ",
         "Install it with: install.packages('rstanarm')",
         call. = FALSE)
  }

  # --- Extract design info from survey.design2 if provided ---
  if (!is.null(design)) {
    design_info <- extract_design(design)
    weights <- design_info$weights
    psu_var <- design_info$cluster
    psu <- as.integer(as.factor(psu_var))
  }

  # --- Auto-detect family ---
  stanreg_family <- stats::family(x)
  family_name <- stanreg_family$family
  family <- switch(family_name,
    "binomial"   = "binomial",
    "gaussian"   = "gaussian",
    stop("Unsupported rstanarm family '", family_name, "'. ",
         "svyder supports: binomial, gaussian.",
         call. = FALSE)
  )

  # --- Extract response variable ---
  y <- x$y
  if (is.null(y)) {
    stop("stanreg object does not contain response data. ",
         "Refit the model to ensure data is stored.", call. = FALSE)
  }

  # --- Extract design matrix ---
  X <- stats::model.matrix(x)
  if (is.null(X)) {
    stop("Could not extract design matrix from stanreg object.",
         call. = FALSE)
  }

  # --- Auto-detect group variable ---
  # rstanarm stores grouping factors in x$glmod
  glmod <- x$glmod
  if (is.null(glmod)) {
    stop("No random effects found in stanreg object. ",
         "der_compute requires a model with group-level effects.",
         call. = FALSE)
  }
  # Get the grouping factor from the first random effect term
  group_list <- glmod$reTrms$flist
  if (length(group_list) == 0L) {
    stop("No grouping factors found in stanreg random effects.",
         call. = FALSE)
  }
  group_factor <- group_list[[1]]
  group <- as.integer(group_factor)

  # --- Auto-detect sigma_theta from posterior ---
  if (is.null(sigma_theta)) {
    draws_mat <- as.matrix(x)
    all_names <- colnames(draws_mat)
    # rstanarm stores group-level SD as Sigma[group:(Intercept),(Intercept)]
    sigma_pars <- all_names[grepl("^Sigma\\[", all_names)]
    if (length(sigma_pars) > 0L) {
      # Sigma parameters in rstanarm are variances, not SDs
      sigma_theta <- sqrt(mean(draws_mat[, sigma_pars[1]]))
    } else {
      stop("Could not auto-detect sigma_theta. ",
           "No 'Sigma[...]' parameters found. ",
           "Please provide sigma_theta manually.",
           call. = FALSE)
    }
  }

  # --- Auto-detect sigma_e for gaussian ---
  if (family == "gaussian" && is.null(sigma_e)) {
    draws_mat <- as.matrix(x)
    all_names <- colnames(draws_mat)
    if ("sigma" %in% all_names) {
      sigma_e <- mean(draws_mat[, "sigma"])
    } else {
      stop("Could not auto-detect sigma_e for gaussian family. ",
           "No 'sigma' parameter found. Please provide sigma_e manually.",
           call. = FALSE)
    }
  }

  # --- Extract draws ---
  extracted <- extract_draws(x)
  draws_mat <- extracted$draws
  param_info <- extracted$param_info

  # Verify dimensions
  p <- ncol(X)
  J <- max(group)

  n_fe <- sum(param_info$type == "fixed")
  n_re <- sum(param_info$type == "random")

  if (n_fe != p) {
    stop(sprintf(
      "Number of fixed effect draws (%d) does not match design matrix columns (%d).",
      n_fe, p
    ), call. = FALSE)
  }

  if (n_re != J) {
    stop(sprintf(
      "Number of random effect draws (%d) does not match number of groups (%d).",
      n_re, J
    ), call. = FALSE)
  }

  # Construct combined draws in svyder order: [beta, theta]
  fe_idx <- which(param_info$type == "fixed")
  re_idx <- which(param_info$type == "random")
  draws_ordered <- draws_mat[, c(fe_idx, re_idx), drop = FALSE]

  # --- Delegate to der_compute.matrix ---
  der_compute.matrix(
    draws_ordered, ...,
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
    design        = NULL  # Already extracted above
  )
}
