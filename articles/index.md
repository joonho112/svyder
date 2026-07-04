# Articles

### Applied track

Workflow-first guides for using svyder in practice — run the pipeline,
choose the variance target, and report results.

- [Getting Started with
  svyder](https://joonho112.github.io/svyder/articles/getting-started.md):

  Run the full compute-classify-correct DER pipeline in about five
  minutes on bundled survey data, then read, tidy, plot, and extract the
  corrected draws.

- [Choosing the Variance
  Target](https://joonho112.github.io/svyder/articles/choosing-the-target.md):

  The one modeling decision the Design Effect Ratio forces you to make:
  what is the sandwich aggregation unit? This vignette contrasts the two
  canonical targets (design-PSU vs model-group) on the bundled
  NSECE-style demo, walks through der_compare(), and explains how
  strata, the finite-cluster correction, the cluster universe, and the
  weight convention all become part of the variance target you declare.

- [The Compute-Classify-Correct
  Pipeline](https://joonho112.github.io/svyder/articles/pipeline.md):

  A step-by-step tour of the svyder workflow: der_compute() builds the
  declared variance target and per-parameter DERs, der_classify() flags
  the parameters that need correction at a threshold tau, and
  der_correct() rescales only those parameters while leaving everything
  else bitwise untouched. Covers der_diagnose() as the one-call
  shortcut, threshold sensitivity, and why selective correction beats a
  blanket adjustment.

- [Case Study: NSECE-Style Survey
  Data](https://joonho112.github.io/svyder/articles/case-study-nsece.md):

  A complete, start-to-finish DER analysis on the bundled synthetic
  NSECE-style survey data: explore the design, run the diagnostic,
  interpret fixed and random effects, read the decomposition, correct
  selectively, compare variance targets, and write up the result for
  publication.

- [Backends: brms, cmdstanr, rstanarm, and
  survey](https://joonho112.github.io/svyder/articles/backends.md):

  How to feed svyder from fitted model objects instead of a raw draws
  matrix. Covers the extract_draws() and extract_design() generics, what
  each backend (brms, cmdstanr, rstanarm, posterior, survey)
  auto-extracts, and what you must always supply yourself — above all,
  the cluster aggregation unit.

### Method track

The mathematical foundations of the DER framework, for readers who want
the derivations.

- [The Design Effect Ratio: Definition and
  Regimes](https://joonho112.github.io/svyder/articles/theory-der.md):

  A methodological introduction to the parameter-level Design Effect
  Ratio (DER): how it generalises the Kish design effect to Bayesian
  hierarchical models, the three regimes it distinguishes, and the
  information-source intuition behind why parameters inherit design
  effects so differently.

- [The Declared Variance Target: Bread and
  Meat](https://joonho112.github.io/svyder/articles/theory-sandwich.md):

  How svyder constructs the design-based sandwich covariance target that
  sits in the numerator of every Design Effect Ratio: the
  penalized-posterior bread, the cluster-aggregated meat, the weight and
  aggregation conventions that make the target a well-defined estimand,
  and the v0.2.0 gaussian score fix.

- [Decomposition Theorems and the Conservation
  Law](https://joonho112.github.io/svyder/articles/theory-decomposition.md):

  The mathematical centerpiece of svyder: under a balanced conjugate
  Gaussian model the Design Effect Ratio admits exact closed forms. This
  article states Theorems 1-2 and the conservation law, and verifies
  them numerically.

- [Classification and Selective
  Correction](https://joonho112.github.io/svyder/articles/theory-classification-correction.md):

  How the Design Effect Ratio turns a per-parameter diagnostic into a
  decision: the four-tier classification of parameters, the threshold
  that separates survey-sensitive parameters from protected ones, and
  the selective block-Cholesky correction (Algorithm 1) that widens only
  the intervals that need widening while leaving shrinkage-protected
  parameters untouched.
