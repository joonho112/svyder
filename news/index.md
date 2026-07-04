# Changelog

## svyder 0.2.0

svyder 0.2.0 makes the **declared variance target** explicit, corrects
the sandwich engine, and ships a full documentation website. The
diagnostic you get now depends on choices you *state* — the aggregation
unit, the weight convention, the strata — rather than on silent
defaults.

### Breaking changes

- [`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md)
  now requires the sandwich aggregation unit explicitly via the
  `cluster` argument (design PSU ids for a design-based target, or
  `group` for a model-aligned target). There is no default: the
  aggregation unit is part of the variance-target definition, and DER
  values can change materially between targets. `psu` is retained as a
  deprecated alias and emits a warning.
- Gaussian score fix: the weighted score residual for the gaussian
  family is now `w * (y - mu) / sigma_e^2`, matching the bread’s working
  weights `w / sigma_e^2`. The previous residual `w * (y - mu)` made the
  sandwich — and hence the DER — scale like `sigma_e^4` under a change
  of response units, so gaussian DERs were wrong by a factor of
  `sigma_e^4` whenever `sigma_e != 1`. Binomial results are unaffected.
- The meat is now stratified, centered within strata, and DF-corrected
  by default (`strata`, `center_meat = TRUE`, `df_correct = TRUE`). Set
  `center_meat = FALSE, df_correct = FALSE` (and `strata = NULL`) to
  reproduce v0.1.0 results exactly.
- [`der_correct()`](https://joonho112.github.io/svyder/reference/der_correct.md)
  methods renamed: the previous per-parameter sqrt(DER) rescaling is now
  `method = "marginal"`; the new default `method = "block_cholesky"`
  transforms the flagged block jointly so its covariance matches the
  sandwich target (the two coincide when exactly one parameter is
  flagged). `method = "cholesky"` was ambiguous and now errors with
  guidance.
- [`der_decompose()`](https://joonho112.github.io/svyder/reference/der_decompose.md)
  computes the protection factor `R_k` structurally from the
  working-weighted within/between split of each covariate; it is no
  longer back-solved from the observed DER, which made the
  predicted-vs-observed comparison circular. Output columns renamed:
  `deff_mean` -\> `deff_used`, `B_mean` -\> `B_used` (now row-specific).
  Requires objects created by svyder \>= 0.2.0 (which store the `$data`
  slots).

### New features

- [`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md)
  gains a `cluster_strata` argument declaring the full cluster universe
  of the design (the stratum of every cluster, indexed by cluster id).
  Clusters with zero sampled observations — e.g., selected PSUs whose
  second-stage sample is empty — then enter the meat as zero score-total
  clusters, contributing to the stratum centering and the `C_h/(C_h-1)`
  correction instead of silently vanishing from the declared target.
  Leave `NULL` when every realized cluster has at least one observation.
- Weight-normalization conventions:
  `normalize = c("unit_mean", "group_size", "none")`. Random-effect DERs
  are not invariant to weight rescaling, so the convention is part of
  the declared target and is recorded in the result.
- Declared-target metadata: results carry `$target` (aggregation unit,
  strata, meat options, weight convention, bread convention), displayed
  by [`print()`](https://rdrr.io/r/base/print.html).
- Tier III exclusion table: results carry `$excluded`, listing
  hyperparameters (`sigma_theta`, and `sigma_e` for gaussian) whose DER
  is undefined, with reasons.
- [`der_compare()`](https://joonho112.github.io/svyder/reference/der_compare.md)
  fast path: given a svyder object (\>= 0.2.0), only the meat is rebuilt
  per target; `y` / `X` / `weights` need not be re-supplied.
- `param_types = NULL` now infers types from the design matrix (columns
  constant within every group are `fe_between`) and reports the
  inference via a message instead of silently defaulting.
- Optional `strata` throughout; singleton strata contribute uncentered
  with a warning.
- Fit-object methods (`brmsfit`, `CmdStanMCMC`, `stanreg`) pass through
  `cluster`, `strata`, `normalize`, `center_meat`, `df_correct`, and
  `beta_prior_sd`; survey designs supply first-stage clusters and
  strata.

### Documentation and website

- New **pkgdown website**: <https://joonho112.github.io/svyder/>.
- A two-track vignette set — **five applied** guides (getting started,
  choosing the variance target, the compute–classify–correct pipeline,
  an NSECE-style case study, and backends) and **four method** articles
  (the DER definition and regimes, the sandwich target, the
  decomposition theorems and conservation law, and classification &
  selective correction).
- roxygen help pages rewritten and expanded: worked, runnable examples
  on the bundled data; the mathematics of the target and the correction
  in `@details`; literature references throughout; and a reference index
  grouped by topic.
- [`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
  [`glance()`](https://generics.r-lib.org/reference/glance.html)
  generics are re-exported, so
  [`library(svyder); tidy(result)`](https://joonho112.github.io/svyder/)
  works without also attaching **generics** (which moves from Suggests
  to Imports).
- New hex logo.

### Bug fixes

- The internal svyder-object validator required a `diagnostics`
  component that no svyder object has ever had; its required-component
  list now matches the actual object structure.

## svyder 0.1.0

- Initial release of the svyder package.
- Core pipeline:
  [`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md),
  [`der_classify()`](https://joonho112.github.io/svyder/reference/der_classify.md),
  [`der_correct()`](https://joonho112.github.io/svyder/reference/der_correct.md),
  [`der_diagnose()`](https://joonho112.github.io/svyder/reference/der_diagnose.md).
- Analysis functions:
  [`der_decompose()`](https://joonho112.github.io/svyder/reference/der_decompose.md),
  [`der_sensitivity()`](https://joonho112.github.io/svyder/reference/der_sensitivity.md),
  [`der_theorem_check()`](https://joonho112.github.io/svyder/reference/der_theorem_check.md),
  [`der_compare()`](https://joonho112.github.io/svyder/reference/der_compare.md).
- Backend support: brms, cmdstanr, rstanarm, survey package.
- Visualization: profile, decomposition, and CI comparison plots.
- Tidy methods:
  [`tidy.svyder()`](https://joonho112.github.io/svyder/reference/tidy.svyder.md),
  [`glance.svyder()`](https://joonho112.github.io/svyder/reference/glance.svyder.md).
- Two bundled datasets: `nsece_demo` and `sim_hlr`.
