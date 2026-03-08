
<!-- README.md is generated from README.Rmd. Please edit that file -->

# svyder <img src="man/figures/logo.png" align="right" height="139" alt="svyder hex sticker" />

<!-- badges: start -->

<!-- badges: end -->

**Design Effect Ratio Diagnostics for Bayesian Survey Models** – *Which
parameters need survey correction?*

## Overview

When Bayesian hierarchical models are fitted to complex survey data, not
all parameters are equally affected by the survey design. Some
parameters (e.g., within-cluster covariates) are highly sensitive to
unequal weighting and clustering, while others (e.g., random effects)
are naturally protected by hierarchical shrinkage.

**svyder** computes parameter-level Design Effect Ratios (DER) that
quantify how much each parameter’s posterior variance should change when
the survey design is properly accounted for. It implements a three-step
compute-classify-correct pipeline: DER values near 1 indicate
well-calibrated posteriors; values substantially above 1 indicate
underestimation of uncertainty.

## Installation

You can install the development version of svyder from
[GitHub](https://github.com/joonho112/svyder):

``` r
# install.packages("pak")
pak::pak("joonho112/svyder")
```

## Quick Demo (5 minutes)

``` r
library(svyder)

# Load bundled NSECE-like survey data
data(nsece_demo)

# Run the full DER diagnostic pipeline
result <- der_diagnose(
  nsece_demo$draws,
  y = nsece_demo$y,
  X = nsece_demo$X,
  group = nsece_demo$group,
  weights = nsece_demo$weights,
  psu = nsece_demo$psu,
  family = "binomial",
  sigma_theta = nsece_demo$sigma_theta,
  param_types = nsece_demo$param_types
)

# Print summary
print(result)
```

The output shows which parameters are flagged for correction:

    svyder diagnostic (54 parameters)
      Family: binomial | N = 6785 | J = 51
      DER range: [0.026, 3.891]
      Threshold (tau): 1.20
      Flagged: 2 / 54 (3.7%)

      Flagged parameters:
        beta[1]              DER = 1.465  [I-b] -> CORRECT
        beta[2]              DER = 3.891  [I-a] -> CORRECT

      Correction applied: 2 parameter(s) rescaled

## DER Profile Plot

Visualize the DER values across all parameters with tier-based coloring:

``` r
plot(result, type = "profile")
```

## Decomposition

Understand *why* each parameter has its observed DER:

``` r
decomp <- der_decompose(result)
head(decomp[, c("param", "param_type", "der", "der_predicted")])
```

## Tidy Output

svyder integrates with the tidyverse via broom-style methods:

``` r
# Per-parameter summary
tidy_df <- tidy(result)
head(tidy_df)

# Model-level summary
glance(result)
```

## Backend Support

svyder works with multiple Bayesian modeling backends:

| Backend  | `der_compute()` method      | Auto-detection                   |
|----------|-----------------------------|----------------------------------|
| Matrix   | `der_compute.matrix()`      | –                                |
| brms     | `der_compute.brmsfit()`     | y, X, group, family, sigma_theta |
| cmdstanr | `der_compute.CmdStanMCMC()` | –                                |
| rstanarm | `der_compute.stanreg()`     | y, X, group, family, sigma_theta |

## Key Functions

| Function              | Purpose                                            |
|-----------------------|----------------------------------------------------|
| `der_diagnose()`      | All-in-one pipeline (compute + classify + correct) |
| `der_compute()`       | Compute DER values                                 |
| `der_classify()`      | Three-tier classification                          |
| `der_correct()`       | Selective Cholesky correction                      |
| `der_decompose()`     | Decompose DER into DEFF, B, R_k, kappa             |
| `der_sensitivity()`   | Threshold sensitivity analysis                     |
| `der_theorem_check()` | Verify theoretical predictions                     |
| `der_compare()`       | Compare across clustering definitions              |

## Citation

If you use svyder in your research, please cite:

> Lee, J. (2026). Design Effect Ratios for Bayesian Survey Models: A
> Diagnostic Framework for Identifying Survey-Sensitive Parameters.
> *arXiv preprint*.

## License

MIT
