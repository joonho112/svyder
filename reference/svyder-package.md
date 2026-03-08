# svyder: Design Effect Ratio Diagnostics for Bayesian Survey Models

Computes parameter-level Design Effect Ratios (DER) for Bayesian
hierarchical models fitted to complex survey data. The package
implements a compute-classify-correct diagnostic framework:

## Details

1.  **Compute**: Calculate DER for each parameter using a sandwich
    variance estimator.

2.  **Classify**: Assign parameters to three tiers based on their
    information source and flag those exceeding a threshold.

3.  **Correct**: Apply selective Cholesky correction to flagged
    parameters only.

## Core pipeline

- [`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md):
  Compute DER values.

- [`der_classify()`](https://joonho112.github.io/svyder/reference/der_classify.md):
  Three-tier classification.

- [`der_correct()`](https://joonho112.github.io/svyder/reference/der_correct.md):
  Selective Cholesky correction.

- [`der_diagnose()`](https://joonho112.github.io/svyder/reference/der_diagnose.md):
  All-in-one wrapper.

## Analysis

- [`der_decompose()`](https://joonho112.github.io/svyder/reference/der_decompose.md):
  Decompose DER into DEFF, B, R_k, kappa.

- [`der_sensitivity()`](https://joonho112.github.io/svyder/reference/der_sensitivity.md):
  Threshold sensitivity analysis.

- [`der_theorem_check()`](https://joonho112.github.io/svyder/reference/der_theorem_check.md):
  Verify theoretical predictions.

- [`der_compare()`](https://joonho112.github.io/svyder/reference/der_compare.md):
  Compare across clustering definitions.

## Backend support

- brms:
  [`der_compute.brmsfit()`](https://joonho112.github.io/svyder/reference/der_compute.md),
  [`extract_draws.brmsfit()`](https://joonho112.github.io/svyder/reference/extract_draws.md)

- cmdstanr:
  [`der_compute.CmdStanMCMC()`](https://joonho112.github.io/svyder/reference/der_compute.md),
  [`extract_draws.CmdStanMCMC()`](https://joonho112.github.io/svyder/reference/extract_draws.md)

- rstanarm:
  [`der_compute.stanreg()`](https://joonho112.github.io/svyder/reference/der_compute.md),
  [`extract_draws.stanreg()`](https://joonho112.github.io/svyder/reference/extract_draws.md)

- survey:
  [`extract_design.survey.design2()`](https://joonho112.github.io/svyder/reference/extract_design.md)

## Bundled datasets

- [nsece_demo](https://joonho112.github.io/svyder/reference/nsece_demo.md):
  Synthetic NSECE-like survey data (binomial, J=51).

- [sim_hlr](https://joonho112.github.io/svyder/reference/sim_hlr.md):
  Balanced hierarchical linear regression (gaussian, J=10).

## See also

Useful links:

- <https://github.com/joonho112/svyder>

- Report bugs at <https://github.com/joonho112/svyder/issues>

## Author

**Maintainer**: JoonHo Lee <jlee296@ua.edu>
([ORCID](https://orcid.org/0009-0006-4019-8703))
