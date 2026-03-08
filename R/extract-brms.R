###############################################################################
# extract-brms.R
# brms integration for the svyder package
# ---------------------------------------------------------------------------
# Provides extract_draws.brmsfit() and der_compute.brmsfit().
# The brms package stores the full model specification (data, formula,
# family, priors) in the fit object, enabling auto-detection of most
# DER computation inputs. Survey weights and PSU must still be provided
# by the user since brms does not store survey design information.
###############################################################################

#' @rdname extract_draws
#'
#' @param pars Character vector of parameter names to extract. If \code{NULL}
#'   (default), extracts all population-level and group-level parameters,
#'   excluding variance/correlation components and diagnostics.
#'
#' @export
extract_draws.brmsfit <- function(x, ..., pars = NULL) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required for extract_draws.brmsfit(). ",
         "Install it with: install.packages('brms')",
         call. = FALSE)
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for extract_draws.brmsfit(). ",
         "Install it with: install.packages('posterior')",
         call. = FALSE)
  }

  # Extract all posterior draws as a draws_matrix
  all_draws <- posterior::as_draws_matrix(x)
  all_names <- colnames(all_draws)

  if (!is.null(pars)) {
    # User-specified parameter selection
    missing_pars <- setdiff(pars, all_names)
    if (length(missing_pars) > 0L) {
      stop("Parameters not found in model: ",
           paste(missing_pars, collapse = ", "),
           call. = FALSE)
    }
    draws_mat <- as.matrix(all_draws[, pars, drop = FALSE])
    param_info <- data.frame(
      original_name = pars,
      svyder_name   = pars,
      type          = ifelse(grepl("^b_", pars), "fixed",
                      ifelse(grepl("^r_", pars), "random", "other")),
      stringsAsFactors = FALSE
    )
  } else {
    # Auto-detect: extract fixed effects (b_) and group-level effects (r_)
    # Exclude: sd_, cor_, sigma, Intercept_sigma, lp__, lprior
    mapped <- .map_brms_params(all_names)

    if (nrow(mapped) == 0L) {
      stop("No population-level (b_) or group-level (r_) parameters found ",
           "in brmsfit object.", call. = FALSE)
    }

    draws_mat <- as.matrix(all_draws[, mapped$original_name, drop = FALSE])
    colnames(draws_mat) <- mapped$svyder_name
    param_info <- mapped
  }

  list(
    draws      = draws_mat,
    param_info = param_info
  )
}


#' Map brms parameter names to svyder convention
#'
#' Converts brms-style parameter names to the svyder naming convention:
#' \itemize{
#'   \item \code{b_Intercept} -> \code{beta[1]}
#'   \item \code{b_x1} -> \code{beta[2]}
#'   \item \code{r_group[level1,Intercept]} -> \code{theta[1]}
#'   \item \code{r_group[level2,Intercept]} -> \code{theta[2]}
#' }
#'
#' Excludes variance components (\code{sd_}, \code{cor_}), residual SD
#' (\code{sigma}), and log-posterior (\code{lp__}, \code{lprior}).
#'
#' @param param_names Character vector of brms parameter names.
#'
#' @return Data frame with columns:
#'   \describe{
#'     \item{original_name}{The original brms parameter name.}
#'     \item{svyder_name}{The mapped svyder parameter name.}
#'     \item{type}{Either "fixed" or "random".}
#'   }
#'
#' @keywords internal
.map_brms_params <- function(param_names) {
  # --- Identify fixed effects (population-level) ---
  fe_names <- param_names[grepl("^b_", param_names)]
  # Exclude distributional parameters (e.g., b_sigma_Intercept)
  fe_names <- fe_names[!grepl("^b_sigma", fe_names)]

  # --- Identify random effects (group-level) ---
  re_names <- param_names[grepl("^r_", param_names)]

  # --- Build mapping ---
  result <- data.frame(
    original_name = character(0),
    svyder_name   = character(0),
    type          = character(0),
    stringsAsFactors = FALSE
  )

  # Map fixed effects: b_Intercept -> beta[1], b_x1 -> beta[2], ...
  if (length(fe_names) > 0L) {
    fe_mapped <- paste0("beta[", seq_along(fe_names), "]")
    result <- rbind(result, data.frame(
      original_name = fe_names,
      svyder_name   = fe_mapped,
      type          = "fixed",
      stringsAsFactors = FALSE
    ))
  }

  # Map random effects: r_group[level,coef] -> theta[1], theta[2], ...
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
#' \strong{brmsfit method}: The brms method auto-detects the response vector
#' \code{y}, design matrix \code{X}, grouping variable \code{group}, model
#' \code{family}, and random effect SD \code{sigma_theta} from the fitted
#' model object. The user must provide \code{weights} (survey weights) and
#' optionally \code{psu} (primary sampling unit indicators), since these are
#' not stored in brms fit objects.
#'
#' The \code{group} argument can optionally be provided to override auto-detection
#' from the random effects structure.
#'
#' @export
der_compute.brmsfit <- function(x, ..., weights, group = NULL,
                                 psu = NULL, family = NULL,
                                 sigma_theta = NULL,
                                 sigma_e = NULL,
                                 beta_prior_sd = 5,
                                 param_types = NULL, design = NULL) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required for der_compute.brmsfit(). ",
         "Install it with: install.packages('brms')",
         call. = FALSE)
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for der_compute.brmsfit(). ",
         "Install it with: install.packages('posterior')",
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
  if (is.null(family)) {
    brms_family <- stats::family(x)
    family_name <- brms_family$family
    # Map brms family names to svyder family names
    family <- switch(family_name,
      "bernoulli"  = "binomial",
      "binomial"   = "binomial",
      "gaussian"   = "gaussian",
      stop("Unsupported brms family '", family_name, "'. ",
           "svyder supports: binomial (bernoulli), gaussian.",
           call. = FALSE)
    )
  }

  # --- Extract model data ---
  model_data <- x$data
  if (is.null(model_data)) {
    stop("brmsfit object does not contain data. ",
         "Refit with 'file_refit = \"on_change\"' or ensure data is stored.",
         call. = FALSE)
  }

  # --- Extract response variable ---
  resp_var <- brms::brmsterms(x$formula)$respform
  resp_name <- all.vars(resp_var)
  if (length(resp_name) == 0L) {
    stop("Could not determine response variable from brms formula.",
         call. = FALSE)
  }
  y <- model_data[[resp_name[1]]]

  # --- Extract design matrix ---
  stan_data <- brms::standata(x)
  X <- stan_data$X
  if (is.null(X)) {
    stop("Could not extract design matrix from brmsfit. ",
         "Ensure the model has population-level (fixed) effects.",
         call. = FALSE)
  }

  # --- Auto-detect group variable ---
  if (is.null(group)) {
    re_terms <- brms::brmsterms(x$formula)$dpars$mu$re
    if (is.null(re_terms) || length(re_terms) == 0L) {
      stop("No random effects found in brmsfit. ",
           "der_compute requires a model with group-level effects.",
           call. = FALSE)
    }
    # Use the first grouping factor
    group_var_name <- re_terms[[1]]$group
    group_factor <- model_data[[group_var_name]]
    if (is.null(group_factor)) {
      stop("Grouping variable '", group_var_name,
           "' not found in model data.", call. = FALSE)
    }
    group <- as.integer(as.factor(group_factor))
  } else {
    group <- as.integer(group)
  }

  # --- Auto-detect sigma_theta from posterior ---
  if (is.null(sigma_theta)) {
    all_draws <- posterior::as_draws_matrix(x)
    all_names <- colnames(all_draws)
    sd_pars <- all_names[grepl("^sd_", all_names)]
    if (length(sd_pars) == 0L) {
      stop("Could not auto-detect sigma_theta. ",
           "No 'sd_' parameters found. Please provide sigma_theta manually.",
           call. = FALSE)
    }
    # Use the first sd parameter (typically sd_group__Intercept)
    sigma_theta <- mean(as.matrix(all_draws[, sd_pars[1]]))
  }

  # --- Auto-detect sigma_e for gaussian ---
  if (family == "gaussian" && is.null(sigma_e)) {
    all_draws <- posterior::as_draws_matrix(x)
    all_names <- colnames(all_draws)
    if ("sigma" %in% all_names) {
      sigma_e <- mean(as.matrix(all_draws[, "sigma"]))
    } else {
      stop("Could not auto-detect sigma_e for gaussian family. ",
           "No 'sigma' parameter found. Please provide sigma_e manually.",
           call. = FALSE)
    }
  }

  # --- Extract draws ---
  extracted <- extract_draws(x)
  draws_mat <- extracted$draws

  # Verify dimensions: should have p fixed + J random columns
  p <- ncol(X)
  J <- max(group)
  d <- p + J

  # Get parameter info
  param_info <- extracted$param_info
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

  # Construct the combined draws matrix in svyder order: [beta, theta]
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
