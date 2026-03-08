###############################################################################
# der_compute.R
# Core entry point for Design Effect Ratio computation
# ---------------------------------------------------------------------------
# S3 generic + matrix method. Other methods (brmsfit, CmdStanMCMC, stanreg)
# will be added in future phases.
###############################################################################

#' Compute Design Effect Ratios
#'
#' Computes Design Effect Ratios (DER) for each parameter in a Bayesian
#' hierarchical model fitted to complex survey data. The DER measures how
#' much each parameter's posterior variance should change when the survey
#' design is properly accounted for.
#'
#' DER is defined as the ratio of the sandwich variance (which accounts for
#' survey design) to the naive MCMC posterior variance. Values near 1 indicate
#' that the posterior is already well-calibrated; values substantially above 1
#' indicate underestimation of uncertainty.
#'
#' @param x A draws matrix (S x d), brmsfit, CmdStanMCMC, or stanreg object.
#'   For the matrix method, columns 1:p are fixed effects and columns
#'   (p+1):(p+J) are random effects.
#' @param ... Additional arguments passed to methods.
#' @param y Response vector (length N).
#' @param X Design matrix (N x p).
#' @param group Integer group indicator (1 to J).
#' @param weights Survey weights (positive, length N).
#' @param psu PSU indicators (default: same as group).
#' @param family Model family: \code{"binomial"} or \code{"gaussian"}.
#' @param sigma_theta Estimated random effect SD.
#' @param sigma_e Residual SD (gaussian only).
#' @param beta_prior_sd Prior SD for fixed effects (default: 5).
#' @param param_types Character vector of length p: \code{"fe_within"} or
#'   \code{"fe_between"}.
#' @param design A \code{survey.design2} object (alternative to weights/psu).
#'
#' @return A \code{svyder} object containing DER values, sandwich matrices,
#'   per-group diagnostics, and original posterior draws.
#'
#' @seealso [der_classify()] for classification, [der_correct()] for
#'   correction, [der_diagnose()] for the all-in-one pipeline.
#' @family core-pipeline
#'
#' @examples
#' data(nsece_demo)
#' result <- der_compute(
#'   nsece_demo$draws,
#'   y = nsece_demo$y, X = nsece_demo$X,
#'   group = nsece_demo$group, weights = nsece_demo$weights,
#'   psu = nsece_demo$psu, family = "binomial",
#'   sigma_theta = nsece_demo$sigma_theta,
#'   param_types = nsece_demo$param_types
#' )
#' print(result)
#'
#' @export
der_compute <- function(x, ...) UseMethod("der_compute")

#' @rdname der_compute
#' @export
der_compute.matrix <- function(x, ..., y, X, group, weights,
                                psu = NULL, family = "binomial",
                                sigma_theta, sigma_e = NULL,
                                beta_prior_sd = 5, param_types = NULL,
                                design = NULL) {

  t0 <- proc.time()
  matched_call <- match.call()

  # --- Extract design info from survey.design2 if provided ---
  if (!is.null(design)) {
    if (!inherits(design, "survey.design2")) {
      stop("'design' must be a survey.design2 object.", call. = FALSE)
    }
    weights <- stats::weights(design)
    psu_var <- design$cluster[[1]]
    psu <- as.integer(as.factor(psu_var))
  }

  # --- Coerce inputs ---
  X     <- as.matrix(X)
  group <- as.integer(group)

  # --- Dimensions ---
  N <- length(y)
  p <- ncol(X)
  J <- max(group)
  d <- p + J

  # --- Validate the draws matrix ---
  if (ncol(x) != d) {
    stop(sprintf(
      "Draws matrix has %d columns but expected %d (p=%d + J=%d).",
      ncol(x), d, p, J
    ), call. = FALSE)
  }

  # --- Split draws ---
  draws_beta  <- x[, seq_len(p), drop = FALSE]
  draws_theta <- x[, (p + 1):d, drop = FALSE]

  # --- Default PSU = group ---
  if (is.null(psu)) {
    psu <- group
  } else {
    psu <- as.integer(psu)
  }

  # --- Default param_types ---
  if (is.null(param_types)) {
    param_types <- rep("fe_between", p)
  }

  # --- Validate inputs ---
  .validate_inputs(y = y, X = X, group = group, weights = weights,
                   family = family, draws_beta = draws_beta,
                   draws_theta = draws_theta)

  # --- Point estimates ---
  beta_hat  <- colMeans(draws_beta)
  theta_hat <- colMeans(draws_theta)

  # --- Linear predictor and fitted values ---
  eta <- as.numeric(X %*% beta_hat) + theta_hat[group]
  mu  <- .compute_mu(family, eta)

  # --- Working weights (with survey weights) ---
  v <- .working_weights(family, mu, weights, sigma_e = sigma_e)

  # --- Weighted residuals ---
  r <- .compute_residuals(family, y, mu, weights)

  # --- Prior precisions ---
  beta_prior_prec  <- if (is.finite(beta_prior_sd)) 1 / beta_prior_sd^2 else 0
  theta_prior_prec <- 1 / sigma_theta^2

  # --- Build sandwich components ---
  H_obs     <- .build_H_obs(X, v, group, J, p, beta_prior_prec, theta_prior_prec)
  H_obs_inv <- .safe_invert(H_obs)
  J_cluster <- .build_J_cluster(X, r, psu, group, p, J)
  V_sand    <- .build_V_sand(H_obs_inv, J_cluster)

  # --- MCMC covariance ---
  draws_all  <- cbind(draws_beta, draws_theta)
  sigma_mcmc <- cov(draws_all)
  sigma_mcmc <- (sigma_mcmc + t(sigma_mcmc)) / 2

  # --- DER = diag(V_sand) / diag(sigma_mcmc) ---
  diag_V    <- diag(V_sand)
  diag_mcmc <- diag(sigma_mcmc)
  stopifnot(all(diag_mcmc > 0))
  der <- diag_V / diag_mcmc

  # --- Parameter names ---
  param_names <- c(paste0("beta[", seq_len(p), "]"),
                   paste0("theta[", seq_len(J), "]"))
  names(der) <- param_names

  # --- Per-group design effect and shrinkage ---
  wt_unweighted <- .working_weights_unweighted(family, mu, sigma_e = sigma_e)
  deff_j <- .compute_deff_j(group, weights)
  B_j    <- .compute_B_j(group, weights, wt_unweighted, sigma_theta^2)

  # --- Placeholder classification and correction fields ---
  all_param_types <- c(param_types, rep("re", J))
  classification <- data.frame(
    param_name = param_names,
    param_type = all_param_types,
    der        = as.numeric(der),
    stringsAsFactors = FALSE
  )

  # --- Timing ---
  elapsed <- (proc.time() - t0)[["elapsed"]]

  # --- Build svyder object ---
  new_svyder(
    der              = der,
    params           = param_names,
    H_obs            = H_obs,
    J_cluster        = J_cluster,
    V_sand           = V_sand,
    sigma_mcmc       = sigma_mcmc,
    deff_j           = deff_j,
    B_j              = B_j,
    classification   = classification,
    tau              = NA_real_,
    corrected_draws  = NULL,
    scale_factors    = rep(1, d),
    original_draws   = draws_all,
    call             = matched_call,
    family           = family,
    n_obs            = N,
    n_groups         = J,
    compute_time     = elapsed
  )
}

#' @rdname der_compute
#' @export
der_compute.default <- function(x, ...) {
  stop("der_compute() does not know how to handle class '",
       class(x)[1], "'.\n",
       "Supported: matrix, brmsfit, CmdStanMCMC, stanreg.",
       call. = FALSE)
}
