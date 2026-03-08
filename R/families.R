###############################################################################
# families.R
# Family-specific computations for the svyder package
# ---------------------------------------------------------------------------
# These are the ONLY place where model-specific math lives.
# All functions are internal (not exported).
###############################################################################

# Dispatcher to family-specific working weight functions
#
# @param family Character string: "binomial" or "gaussian"
# @param mu Numeric vector of fitted means (length N)
# @param w Numeric vector of survey weights (length N)
# @param sigma_e Numeric scalar, residual SD (required for gaussian)
# @return Numeric vector of working weights (length N)
.working_weights <- function(family, mu, w, sigma_e = NULL) {
  switch(family,
    binomial = .working_weights_binomial(mu, w),
    gaussian = .working_weights_gaussian(w, sigma_e),
    stop(sprintf("Unknown family '%s'. Supported families: 'binomial', 'gaussian'.", family))
  )
}

# Binomial working weights: v_i = w_i * mu_i * (1 - mu_i)
# Ported from standalone compute_der.R lines 62-63
#
# @param mu Numeric vector of fitted probabilities (length N)
# @param w Numeric vector of survey weights (length N)
# @return Numeric vector of working weights (length N)
.working_weights_binomial <- function(mu, w) {
  w * mu * (1 - mu)
}

# Gaussian working weights: v_i = w_i / sigma_e^2
#
# @param w Numeric vector of survey weights (length N)
# @param sigma_e Numeric scalar, residual standard deviation
# @return Numeric vector of working weights (length N)
.working_weights_gaussian <- function(w, sigma_e) {
  if (is.null(sigma_e) || !is.numeric(sigma_e) || length(sigma_e) != 1L || sigma_e <= 0) {
    stop("'sigma_e' must be a positive numeric scalar for the gaussian family.")
  }

  w / sigma_e^2
}

# Link inverse function: compute mu from linear predictor eta
#
# @param family Character string: "binomial" or "gaussian"
# @param eta Numeric vector of linear predictor values
# @return Numeric vector of fitted means (mu)
.compute_mu <- function(family, eta) {
  switch(family,
    binomial = 1 / (1 + exp(-eta)),
    gaussian = eta,
    stop(sprintf("Unknown family '%s'. Supported families: 'binomial', 'gaussian'.", family))
  )
}

# Weighted residuals for score computation
# Ported from standalone compute_der.R line 98
#
# For both families: r_i = w_i * (y_i - mu_i)
#
# @param family Character string: "binomial" or "gaussian"
# @param y Numeric vector of responses (length N)
# @param mu Numeric vector of fitted means (length N)
# @param w Numeric vector of survey weights (length N)
# @return Numeric vector of weighted residuals (length N)
.compute_residuals <- function(family, y, mu, w) {
  # Same formula for both families
  w * (y - mu)
}

# Working weights WITHOUT survey weights (for shrinkage computation)
# Returns mu*(1-mu) for binomial, 1/sigma_e^2 for gaussian
#
# @param family Character string: "binomial" or "gaussian"
# @param mu Numeric vector of fitted means (length N)
# @param sigma_e Numeric scalar, residual SD (required for gaussian)
# @return Numeric vector of unweighted working weights (length N)
.working_weights_unweighted <- function(family, mu, sigma_e = NULL) {
  switch(family,
    binomial = mu * (1 - mu),
    gaussian = {
      if (is.null(sigma_e) || !is.numeric(sigma_e) || length(sigma_e) != 1L || sigma_e <= 0) {
        stop("'sigma_e' must be a positive numeric scalar for the gaussian family.")
      }
      rep(1 / sigma_e^2, length(mu))
    },
    stop(sprintf("Unknown family '%s'. Supported families: 'binomial', 'gaussian'.", family))
  )
}
