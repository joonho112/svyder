# Package index

## Core pipeline

The compute-classify-correct workflow. Use
[`der_diagnose()`](https://joonho112.github.io/svyder/reference/der_diagnose.md)
for the all-in-one call, or invoke the three steps individually for
finer control.

- [`der_diagnose()`](https://joonho112.github.io/svyder/reference/der_diagnose.md)
  : All-in-One DER Diagnostic
- [`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md)
  : Compute Design Effect Ratios
- [`der_classify()`](https://joonho112.github.io/svyder/reference/der_classify.md)
  : Classify Parameters by Design Sensitivity
- [`der_correct()`](https://joonho112.github.io/svyder/reference/der_correct.md)
  : Apply Selective Correction to Flagged Parameters

## Analysis

Understand, verify, and re-target DER results.

- [`der_decompose()`](https://joonho112.github.io/svyder/reference/der_decompose.md)
  : Decompose DER into Components
- [`der_sensitivity()`](https://joonho112.github.io/svyder/reference/der_sensitivity.md)
  : Sensitivity Analysis Across Threshold Values
- [`der_theorem_check()`](https://joonho112.github.io/svyder/reference/der_theorem_check.md)
  : Verify Theoretical DER Predictions
- [`der_compare()`](https://joonho112.github.io/svyder/reference/der_compare.md)
  : Compare DER Across Variance Targets

## Object methods

S3 methods for `svyder` objects: printing, summaries, corrected draws,
and broom-style tidiers.

- [`print(`*`<svyder>`*`)`](https://joonho112.github.io/svyder/reference/print.svyder.md)
  : Print a svyder Object
- [`summary(`*`<svyder>`*`)`](https://joonho112.github.io/svyder/reference/summary.svyder.md)
  : Summarize a svyder Object
- [`tidy(`*`<svyder>`*`)`](https://joonho112.github.io/svyder/reference/tidy.svyder.md)
  : Tidy a svyder Object
- [`glance(`*`<svyder>`*`)`](https://joonho112.github.io/svyder/reference/glance.svyder.md)
  : Glance at a svyder Object
- [`as.matrix(`*`<svyder>`*`)`](https://joonho112.github.io/svyder/reference/as.matrix.svyder.md)
  : Extract Draws Matrix from a svyder Object
- [`is.svyder()`](https://joonho112.github.io/svyder/reference/is.svyder.md)
  : Test if an Object is a svyder Object

## Visualization

Diagnostic plots for DER results: profile, decomposition, and comparison
views, with an automatic ggplot2 / base R fallback.

- [`plot(`*`<svyder>`*`)`](https://joonho112.github.io/svyder/reference/plot.svyder.md)
  : Plot DER Diagnostic Results
- [`autoplot.svyder()`](https://joonho112.github.io/svyder/reference/autoplot.svyder.md)
  : Create a ggplot2 Visualization of DER Results

## Backends

Generic extractors for posterior draws and survey design objects.
Methods exist for brms, cmdstanr, rstanarm, posterior, and the survey
package.

- [`extract_draws()`](https://joonho112.github.io/svyder/reference/extract_draws.md)
  : Extract Posterior Draws from Model Objects
- [`extract_design()`](https://joonho112.github.io/svyder/reference/extract_design.md)
  : Extract Survey Design Information

## Bundled datasets

Synthetic datasets shipping pre-computed posterior draws, so no Stan
installation is required.

- [`nsece_demo`](https://joonho112.github.io/svyder/reference/nsece_demo.md)
  : Synthetic NSECE-Like Survey Data
- [`sim_hlr`](https://joonho112.github.io/svyder/reference/sim_hlr.md) :
  Simulated Hierarchical Linear Regression Data

## Package

- [`svyder`](https://joonho112.github.io/svyder/reference/svyder-package.md)
  [`svyder-package`](https://joonho112.github.io/svyder/reference/svyder-package.md)
  : svyder: Design Effect Ratio Diagnostics for Bayesian Survey Models
