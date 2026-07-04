###############################################################################
# der_compute.R
# Core entry point for Design Effect Ratio computation
# ---------------------------------------------------------------------------
# S3 generic + matrix method. Fit-object methods (brmsfit, CmdStanMCMC,
# stanreg) live in extract-*.R and delegate here.
###############################################################################

#' Compute Design Effect Ratios
#'
#' Computes Design Effect Ratios (DER) for each parameter in a Bayesian
#' hierarchical model fitted to complex survey data. The DER of a parameter
#' is the ratio of its design-based sandwich variance -- computed for a
#' \emph{declared variance target} -- to its model-based MCMC posterior
#' variance. Values substantially above 1 indicate that the posterior
#' understates design-based uncertainty; values at or below 1 indicate that
#' correction is unnecessary (and typically harmful).
#'
#' @section The declared variance target:
#' The sandwich variance is not unique: it depends on choices that are part
#' of the estimand, not implementation details. \code{der_compute()} therefore
#' requires them explicitly and records them in the returned object:
#' \itemize{
#'   \item \code{cluster} -- the aggregation unit of the score outer product
#'     (the design PSU for a design-based target, or the model group for a
#'     model-aligned target). There is no default: DER values can change
#'     materially between the two (in the NSECE application the flagged set
#'     moves from 1 of 54 to 28 of 54, including 27 of the 51 state
#'     effects).
#'   \item \code{strata} -- optional design strata; cluster totals are
#'     centered within strata and a \eqn{C_h/(C_h-1)} correction is applied.
#'   \item \code{normalize} -- the survey-weight convention. Random-effect
#'     DERs are not invariant to weight rescaling, so the convention is part
#'     of the target. \code{"unit_mean"} (default) rescales weights to mean 1
#'     overall; \code{"group_size"} rescales within each group to sum to the
#'     group size; \code{"none"} uses the weights as supplied.
#' }
#'
#' @section The bread:
#' The bread is the penalized (posterior) Hessian at the posterior mean: the
#' pseudo-likelihood curvature plus the random-effect prior precision
#' \eqn{1/\sigma_\theta^2}. The likelihood-only Hessian is singular for any
#' model with an intercept and a full set of group effects, so the prior term
#' is what makes the sandwich well-defined. By the diffuse-beta convention
#' the fixed-effect prior contributes no curvature by default
#' (\code{beta_prior_sd = Inf}); supply a finite value to reproduce pipelines
#' that used an informative beta prior in the bread.
#'
#' @param x A draws matrix (S x d), brmsfit, CmdStanMCMC, or stanreg object.
#'   For the matrix method, columns 1:p are fixed effects and columns
#'   (p+1):(p+J) are random effects.
#' @param ... Additional arguments passed to methods.
#' @param y Response vector (length N).
#' @param X Design matrix (N x p).
#' @param group Integer group indicator (1 to J).
#' @param weights Survey weights (positive, length N).
#' @param cluster Sandwich aggregation unit (length N). \strong{Required}
#'   (unless supplied via \code{design}): pass the design PSU ids for a
#'   design-based target or \code{group} for a model-aligned target.
#' @param strata Optional stratum indicator (length N).
#' @param cluster_strata Optional vector giving the stratum of \emph{every}
#'   cluster in the declared design universe, indexed by cluster id (length
#'   G, the universe size). Supplying it declares the cluster universe
#'   explicitly, so cluster ids with no sampled observations -- selected
#'   PSUs whose second-stage (e.g., Poisson) sample is empty -- enter the
#'   meat as zero score-total clusters, contributing to the stratum
#'   centering and the \eqn{C_h/(C_h-1)} correction rather than silently
#'   vanishing from the target. When supplied, \code{cluster} must contain
#'   integer ids in \code{1:length(cluster_strata)} and \code{strata} is
#'   required. Leave \code{NULL} when every cluster of the realized design
#'   has at least one observation (e.g., the NSECE application).
#' @param family Model family: \code{"binomial"} or \code{"gaussian"}.
#' @param sigma_theta Estimated random effect SD (plug-in).
#' @param sigma_e Residual SD (gaussian only).
#' @param normalize Weight convention: \code{"unit_mean"} (default),
#'   \code{"group_size"}, or \code{"none"}. Recorded in the result.
#' @param center_meat Center cluster score totals within strata
#'   (default \code{TRUE}). Set \code{FALSE} only to reproduce v1 results.
#' @param df_correct Apply the \eqn{C_h/(C_h-1)} stratum correction
#'   (default \code{TRUE}). Set \code{FALSE} only to reproduce v1 results.
#' @param beta_prior_sd Prior SD for fixed effects in the bread
#'   (default \code{Inf} = diffuse-beta convention, no ridge).
#' @param param_types Character vector of length p: \code{"fe_within"} or
#'   \code{"fe_between"}. If \code{NULL}, types are inferred from the design
#'   matrix (columns constant within every group are \code{"fe_between"})
#'   and the inference is reported via a message.
#' @param design A \code{survey.design2} object; supplies \code{weights},
#'   \code{cluster} (first-stage), and \code{strata} (first-stage).
#' @param psu Deprecated alias for \code{cluster}.
#'
#' @return A \code{svyder} object (a structured list) with, among others:
#'   \describe{
#'     \item{der}{Named numeric vector of Design Effect Ratios, one per
#'       parameter (\code{beta[1..p]} then \code{theta[1..J]}).}
#'     \item{V_sand}{Sandwich (declared-target) variance matrix
#'       \eqn{H^{-1} J_c H^{-1}} (\eqn{d \times d}).}
#'     \item{sigma_mcmc}{Model-based MCMC posterior covariance
#'       (\eqn{d \times d}); \code{der} is its diagonal ratio to
#'       \code{V_sand}.}
#'     \item{target}{List recording the declared variance target: number of
#'       clusters, strata, centering / DF-correction flags, weight
#'       convention, and bread convention.}
#'     \item{excluded}{Data frame of Tier III parameters (hyperparameters)
#'       excluded from the DER diagnostic, with reasons.}
#'     \item{data}{List of data slots (\code{X}, \code{group}, \code{weights},
#'       \code{cluster}, \code{strata}, working weights \code{v}, score
#'       residuals \code{r}) retained for re-targeting and decomposition.}
#'     \item{classification}{Data frame with one row per parameter
#'       (\code{param_name}, \code{param_type}, \code{der}); tiers are filled
#'       in by [der_classify()].}
#'     \item{deff_j, B_j}{Per-group Kish design effects and shrinkage factors
#'       (length J).}
#'   }
#'
#' @references
#' Lee, J., Williams, M. R., & Savitsky, T. D. (2026). Design Effect Ratios
#' for Bayesian Survey Models: A Diagnostic Framework for Identifying
#' Survey-Sensitive Parameters. \emph{Journal of Survey Statistics and
#' Methodology}. Submitted.
#'
#' Savitsky, T. D., & Williams, M. R. (2022). Pseudo Bayesian mixed models
#' under informative sampling. \emph{Journal of Official Statistics}, 38(4),
#' 901--928. \doi{10.2478/jos-2022-0039}
#'
#' Williams, M. R., & Savitsky, T. D. (2021). Uncertainty estimation for
#' pseudo-Bayesian inference under complex sampling. \emph{International
#' Statistical Review}, 89(1), 72--107. \doi{10.1111/insr.12376}
#'
#' Binder, D. A. (1983). On the variances of asymptotically normal estimators
#' from complex surveys. \emph{International Statistical Review}, 51(3),
#' 279--292. \doi{10.2307/1402588}
#'
#' Kish, L. (1965). \emph{Survey Sampling}. John Wiley & Sons.
#'
#' @seealso [der_classify()] for classification, [der_correct()] for
#'   correction, [der_compare()] for comparing targets,
#'   [der_diagnose()] for the all-in-one pipeline.
#' @family core-pipeline
#'
#' @examples
#' data(nsece_demo)
#' result <- der_compute(
#'   nsece_demo$draws,
#'   y = nsece_demo$y, X = nsece_demo$X,
#'   group = nsece_demo$group, weights = nsece_demo$weights,
#'   cluster = nsece_demo$psu, family = "binomial",
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
                               cluster, strata = NULL,
                               cluster_strata = NULL,
                               family = "binomial",
                               sigma_theta, sigma_e = NULL,
                               normalize = c("unit_mean", "group_size", "none"),
                               center_meat = TRUE, df_correct = TRUE,
                               beta_prior_sd = Inf, param_types = NULL,
                               design = NULL, psu = NULL) {

  t0 <- proc.time()
  matched_call <- match.call()
  normalize <- match.arg(normalize)

  # --- Deprecated psu alias ---
  if (!is.null(psu)) {
    if (missing(cluster)) {
      warning("'psu' is deprecated; use 'cluster'. ",
              "Treating psu as the cluster (aggregation unit).",
              call. = FALSE)
      cluster <- psu
    } else {
      stop("Supply either 'cluster' or the deprecated 'psu', not both.",
           call. = FALSE)
    }
  }

  # --- Extract design info from survey.design2 if provided ---
  if (!is.null(design)) {
    if (!inherits(design, "survey.design2")) {
      stop("'design' must be a survey.design2 object.", call. = FALSE)
    }
    weights <- stats::weights(design)
    cluster_var <- design$cluster[[1]]
    cluster <- as.integer(as.factor(cluster_var))
    if (is.null(strata) && !is.null(design$strata)) {
      strata_var <- design$strata[[1]]
      if (length(unique(strata_var)) > 1L) {
        strata <- as.integer(as.factor(strata_var))
      }
    }
  }

  # --- The aggregation unit is part of the target definition: no default ---
  if (missing(cluster) || is.null(cluster)) {
    stop("'cluster' must be supplied explicitly -- the sandwich aggregation ",
         "unit is part of the variance-target definition.\n",
         "  * design-based target: cluster = <design PSU ids>\n",
         "  * model-aligned target: cluster = group\n",
         "See der_compare() to report both.", call. = FALSE)
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

  # --- Cluster / strata coercion and validation ---
  if (length(cluster) != N) {
    stop(sprintf(
      "Dimension mismatch: length(y) = %d but length(cluster) = %d.",
      N, length(cluster)
    ), call. = FALSE)
  }
  if (anyNA(cluster)) {
    stop("'cluster' contains missing values.", call. = FALSE)
  }
  if (!is.null(cluster_strata)) {
    # Declared cluster universe: 'cluster' must already be integer ids into
    # 1:length(cluster_strata); ids absent from the data (selected-but-empty
    # clusters, e.g. PSUs whose second-stage Poisson sample is empty) are
    # legitimate and contribute zero score totals to the meat.
    cluster_int <- as.integer(cluster)
    G <- length(cluster_strata)
    if (anyNA(cluster_int) || any(cluster_int < 1L) || any(cluster_int > G)) {
      stop("With 'cluster_strata', 'cluster' must contain integer ids in ",
           "1:length(cluster_strata).", call. = FALSE)
    }
  } else {
    cluster_int <- as.integer(as.factor(cluster))
    G <- max(cluster_int)
  }

  if (!is.null(strata)) {
    if (length(strata) != N) {
      stop(sprintf(
        "Dimension mismatch: length(y) = %d but length(strata) = %d.",
        N, length(strata)
      ), call. = FALSE)
    }
    if (anyNA(strata)) {
      stop("'strata' contains missing values.", call. = FALSE)
    }
    strata_int <- as.integer(as.factor(strata))
    H_n <- max(strata_int)
  } else {
    strata_int <- NULL
    H_n <- 1L
  }

  cs_int <- NULL
  if (!is.null(cluster_strata)) {
    if (is.null(strata)) {
      stop("'cluster_strata' requires 'strata'.", call. = FALSE)
    }
    cs_int <- match(cluster_strata, sort(unique(strata)))
    if (anyNA(cs_int)) {
      stop("Every value of 'cluster_strata' must appear in 'strata' ",
           "(a stratum consisting entirely of empty clusters is not ",
           "estimable).", call. = FALSE)
    }
    bad <- strata_int != cs_int[cluster_int]
    if (any(bad)) {
      stop(sprintf(
        "'strata' disagrees with 'cluster_strata' at %d observation(s).",
        sum(bad)), call. = FALSE)
    }
  }

  # --- param_types: explicit or inferred (never a silent default) ---
  if (is.null(param_types)) {
    param_types <- .infer_param_types(X, group)
    message(sprintf(
      "param_types inferred from the design matrix: %s. Supply 'param_types' to override.",
      paste(sprintf("%s=%s", colnames(X) %||% paste0("X", seq_len(p)),
                    param_types), collapse = ", ")
    ))
  }
  if (length(param_types) != p ||
      !all(param_types %in% c("fe_within", "fe_between"))) {
    stop("'param_types' must be length ncol(X) with values 'fe_within' or 'fe_between'.",
         call. = FALSE)
  }

  # --- Validate inputs ---
  .validate_inputs(y = y, X = X, group = group, weights = weights,
                   family = family, draws_beta = draws_beta,
                   draws_theta = draws_theta)

  # --- Weight convention (part of the target; recorded) ---
  weights_raw_mean <- mean(weights)
  if (normalize == "unit_mean") {
    weights <- weights / mean(weights)
  } else if (normalize == "group_size") {
    for (j in seq_len(J)) {
      idx_j <- which(group == j)
      weights[idx_j] <- weights[idx_j] * length(idx_j) / sum(weights[idx_j])
    }
  }

  # --- Point estimates (posterior-mean plug-in) ---
  beta_hat  <- colMeans(draws_beta)
  theta_hat <- colMeans(draws_theta)

  # --- Linear predictor and fitted values ---
  eta <- as.numeric(X %*% beta_hat) + theta_hat[group]
  mu  <- .compute_mu(family, eta)

  # --- Working weights (with survey weights) ---
  v <- .working_weights(family, mu, weights, sigma_e = sigma_e)

  # --- Weighted score residuals ---
  r <- .compute_residuals(family, y, mu, weights, sigma_e = sigma_e)

  # --- Prior precisions (bread convention) ---
  beta_prior_prec  <- if (is.finite(beta_prior_sd)) 1 / beta_prior_sd^2 else 0
  theta_prior_prec <- 1 / sigma_theta^2

  # --- Build sandwich components ---
  H_obs     <- .build_H_obs(X, v, group, J, p, beta_prior_prec, theta_prior_prec)
  H_obs_inv <- .safe_invert(H_obs)
  J_cluster <- .build_J_cluster(X, r, cluster_int, group, p, J,
                                strata = strata_int,
                                cluster_strata = cs_int,
                                center = center_meat,
                                df_correct = df_correct)
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

  # --- Placeholder classification (tiers assigned by der_classify) ---
  all_param_types <- c(param_types, rep("re", J))
  classification <- data.frame(
    param_name = param_names,
    param_type = all_param_types,
    der        = as.numeric(der),
    stringsAsFactors = FALSE
  )

  # --- Tier III: hyperparameters are outside the declared DER target ---
  excluded_params <- "sigma_theta"
  excluded_reasons <- "random-effect prior hyperparameter outside the data-level phi=(beta, theta) score block; DER undefined"
  if (family == "gaussian") {
    excluded_params <- c(excluded_params, "sigma_e")
    excluded_reasons <- c(
      excluded_reasons,
      "gaussian dispersion plug-in; svyder conditions on sigma_e and constructs the location-block sandwich for phi=(beta, theta); DER undefined"
    )
  }
  excluded <- data.frame(
    param  = excluded_params,
    tier   = rep("III", length(excluded_params)),
    reason = excluded_reasons,
    stringsAsFactors = FALSE
  )

  # --- Declared target metadata ---
  target <- list(
    cluster_n        = G,
    cluster_n_empty  = if (!is.null(cs_int))
      G - length(unique(cluster_int)) else 0L,
    cluster_universe_declared = !is.null(cs_int),
    cluster_is_group = isTRUE(all.equal(as.integer(cluster_int),
                                        as.integer(as.factor(group)))),
    strata_n         = H_n,
    center           = center_meat,
    df_correct       = df_correct,
    normalize        = normalize,
    weights_raw_mean = weights_raw_mean,
    plug_in          = "posterior_mean",
    beta_prior_sd    = beta_prior_sd,
    sigma_theta      = sigma_theta
  )

  # --- Data slots for re-targeting and decomposition ---
  data_slots <- list(
    X       = X,
    group   = group,
    weights = weights,
    cluster = cluster_int,
    strata  = strata_int,
    cluster_strata = cs_int,
    v       = v,
    r       = r
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
    compute_time     = elapsed,
    target           = target,
    excluded         = excluded,
    data             = data_slots
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


# Infer fixed-effect parameter types from the design matrix
#
# A column that is constant within every group carries only between-group
# identifying variation ("fe_between"); any within-group variation makes it
# "fe_within". The intercept is constant within groups, hence "fe_between".
#
# @param X Design matrix (N x p)
# @param group Integer group indicator (length N)
# @return Character vector of length p
.infer_param_types <- function(X, group) {
  p <- ncol(X)
  types <- character(p)
  for (k in seq_len(p)) {
    within_var <- tapply(X[, k], group, function(z) max(z) - min(z))
    types[k] <- if (all(within_var < 1e-12)) "fe_between" else "fe_within"
  }
  types
}

# Null-coalescing helper
`%||%` <- function(a, b) if (is.null(a)) b else a
