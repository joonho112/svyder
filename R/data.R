###############################################################################
# data.R
# Dataset documentation for bundled svyder datasets
###############################################################################

#' Synthetic NSECE-Like Survey Data
#'
#' A synthetic dataset mimicking the NSECE 2019 survey structure
#' for demonstration of the DER diagnostic pipeline. Contains
#' N = 6785 observations across J = 51 states with unequal survey
#' weights, clustered PSU structure, and three fixed-effect covariates
#' (intercept, within-cluster poverty, between-cluster tiered
#' reimbursement policy).
#'
#' @format A list with components:
#' \describe{
#'   \item{draws}{Matrix of posterior draws (4000 x 54), columns 1:3 are
#'     fixed effects (beta), columns 4:54 are random effects (theta).}
#'   \item{y}{Binary outcome vector (length 6785).}
#'   \item{X}{Design matrix (6785 x 3) with columns: intercept,
#'     poverty_cwc (group-mean centered), tiered_reim (binary policy).}
#'   \item{group}{Integer state group indicator (1 to 51).}
#'   \item{weights}{Survey weights (positive, length 6785). Log-normal
#'     distributed, normalized within state.}
#'   \item{psu}{PSU indicators (integer, length 6785).}
#'   \item{param_types}{Character vector of length 3:
#'     \code{c("fe_between", "fe_within", "fe_between")}.}
#'   \item{family}{Model family: \code{"binomial"}.}
#'   \item{sigma_theta}{Random effect SD (0.66).}
#'   \item{N}{Number of observations (6785).}
#'   \item{J}{Number of groups (51).}
#'   \item{p}{Number of fixed effects (3).}
#' }
#'
#' @source Synthetic data generated to mimic NSECE 2019 structure.
#'   See \code{data-raw/generate_nsece_demo.R}.
#'
#' @examples
#' data(nsece_demo)
#' str(nsece_demo, max.level = 1)
"nsece_demo"


#' Simulated Hierarchical Linear Regression Data
#'
#' A small balanced Gaussian hierarchical model dataset for quick
#' testing and demonstration. Contains J = 10 groups with n_j = 20
#' observations each (N = 200 total), equal weights (DEFF = 1),
#' and two fixed-effect covariates (intercept + within-cluster
#' covariate).
#'
#' With equal weights the design effect is 1.0, so DER values should
#' be close to 1.0 across all parameters. This dataset is useful for
#' verifying that the pipeline correctly identifies the absence of
#' design effects.
#'
#' @format A list with components:
#' \describe{
#'   \item{draws}{Matrix of posterior draws (4000 x 12), columns 1:2 are
#'     fixed effects (beta), columns 3:12 are random effects (theta).}
#'   \item{y}{Continuous outcome vector (length 200).}
#'   \item{X}{Design matrix (200 x 2) with columns: intercept, x_within.}
#'   \item{group}{Integer group indicator (1 to 10).}
#'   \item{weights}{Survey weights (all 1.0, length 200).}
#'   \item{psu}{PSU indicators (same as group).}
#'   \item{param_types}{Character vector of length 2:
#'     \code{c("fe_between", "fe_within")}.}
#'   \item{family}{Model family: \code{"gaussian"}.}
#'   \item{sigma_theta}{Random effect SD (0.5).}
#'   \item{sigma_e}{Residual SD (1.0).}
#'   \item{N}{Number of observations (200).}
#'   \item{J}{Number of groups (10).}
#'   \item{p}{Number of fixed effects (2).}
#'   \item{B_ref}{Analytical shrinkage factor (5/6).}
#'   \item{deff_ref}{Reference design effect (1.0).}
#' }
#'
#' @source Synthetic data. See \code{data-raw/generate_sim_hlr.R}.
#'
#' @examples
#' data(sim_hlr)
#' str(sim_hlr, max.level = 1)
"sim_hlr"
