# Articles

### Getting Started

- [Getting Started with
  svyder](https://joonho112.github.io/svyder/articles/getting-started.md):

  A hands-on introduction to Design Effect Ratio diagnostics with
  svyder: run the full pipeline in one call, interpret the three-tier
  classification, visualize DER profiles, extract corrected posterior
  draws, and build a pipe-friendly workflow — all using bundled
  datasets.

### Theory

- [Understanding Design Effect
  Ratios](https://joonho112.github.io/svyder/articles/understanding-der.md):

  A deep dive into the Design Effect Ratio (DER) framework: definition
  and three regimes, why parameters differ in design sensitivity, the
  decomposition theorems linking DER to Kish DEFF and hierarchical
  shrinkage, the conservation law, three-tier classification, threshold
  sensitivity, and comparison between complex survey and equal-weight
  data.

- [Decomposition Theorems: Why DER Differs Across
  Parameters](https://joonho112.github.io/svyder/articles/decomposition-theorems.md):

  A deep dive into the DER decomposition theorems: understand why
  different parameters have different design sensitivities, verify the
  closed-form approximations against empirical values, and explore the
  conservation law that links hierarchical shrinkage to survey design
  effects.

### Applications

- [The Compute-Classify-Correct
  Pipeline](https://joonho112.github.io/svyder/articles/pipeline-in-depth.md):

  A comprehensive walkthrough of the three-step DER diagnostic pipeline:
  compute sandwich-based DER values, classify parameters into tiers with
  threshold-based flagging, and apply selective Cholesky correction to
  flagged parameters only. Includes sensitivity analysis,
  cross-clustering comparison, and practical guidance for working with
  corrected posterior draws.

- [Case Study: NSECE Survey Data
  Analysis](https://joonho112.github.io/svyder/articles/case-study-nsece.md):

  A comprehensive case study demonstrating the full DER diagnostic
  workflow on NSECE-like survey data: explore the survey design, run
  diagnostics, interpret fixed and random effect results, verify
  decomposition theorems, compare selective vs blanket correction, and
  report results for publication.
