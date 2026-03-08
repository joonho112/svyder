#' svyder: Design Effect Ratio Diagnostics for Bayesian Survey Models
#'
#' Computes parameter-level Design Effect Ratios (DER) for Bayesian
#' hierarchical models fitted to complex survey data. The package implements
#' a compute-classify-correct diagnostic framework:
#'
#' \enumerate{
#'   \item \strong{Compute}: Calculate DER for each parameter using a sandwich
#'     variance estimator.
#'   \item \strong{Classify}: Assign parameters to three tiers based on their
#'     information source and flag those exceeding a threshold.
#'   \item \strong{Correct}: Apply selective Cholesky correction to flagged
#'     parameters only.
#' }
#'
#' @section Core pipeline:
#' \itemize{
#'   \item [der_compute()]: Compute DER values.
#'   \item [der_classify()]: Three-tier classification.
#'   \item [der_correct()]: Selective Cholesky correction.
#'   \item [der_diagnose()]: All-in-one wrapper.
#' }
#'
#' @section Analysis:
#' \itemize{
#'   \item [der_decompose()]: Decompose DER into DEFF, B, R_k, kappa.
#'   \item [der_sensitivity()]: Threshold sensitivity analysis.
#'   \item [der_theorem_check()]: Verify theoretical predictions.
#'   \item [der_compare()]: Compare across clustering definitions.
#' }
#'
#' @section Backend support:
#' \itemize{
#'   \item brms: [der_compute.brmsfit()], [extract_draws.brmsfit()]
#'   \item cmdstanr: [der_compute.CmdStanMCMC()], [extract_draws.CmdStanMCMC()]
#'   \item rstanarm: [der_compute.stanreg()], [extract_draws.stanreg()]
#'   \item survey: [extract_design.survey.design2()]
#' }
#'
#' @section Bundled datasets:
#' \itemize{
#'   \item [nsece_demo]: Synthetic NSECE-like survey data (binomial, J=51).
#'   \item [sim_hlr]: Balanced hierarchical linear regression (gaussian, J=10).
#' }
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom stats cov family model.matrix quantile qnorm sd var weights
NULL

# For ggplot2 .data pronoun (used in aes())
utils::globalVariables(".data")
