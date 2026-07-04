# The Compute-Classify-Correct Pipeline

## Overview

svyder implements a three-step workflow for correcting Bayesian survey
models: **compute** the parameter-level Design Effect Ratios (DER),
**classify** each parameter into a correction tier, and **correct** only
the parameters that need it. This vignette is for practitioners who want
to understand what happens at each step and when to run the steps
separately instead of all at once. You will see how the correction is
*selective by construction* — flagged parameters are rescaled and
everything else is left exactly as it was.

The single call
[`der_diagnose()`](https://joonho112.github.io/svyder/reference/der_diagnose.md)
runs all three steps for you. But separating them gives you control: you
can inspect the DERs before choosing a threshold, re-run classification
at several thresholds, or correct with an alternative method. We work
with the built-in `nsece_demo` dataset throughout, so every chunk runs
with no Stan installation.

We cover, in order:

- **Step 1 — Compute.**
  [`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md):
  what it needs and what it returns.
- **Step 2 — Classify.**
  [`der_classify()`](https://joonho112.github.io/svyder/reference/der_classify.md):
  the tiers and the threshold `tau`.
- **Step 3 — Correct.**
  [`der_correct()`](https://joonho112.github.io/svyder/reference/der_correct.md):
  block-Cholesky vs marginal, and the bitwise-identity guarantee for
  unflagged parameters.
- **All in one.**
  [`der_diagnose()`](https://joonho112.github.io/svyder/reference/der_diagnose.md)
  as the equivalent shortcut.
- **Threshold sensitivity.**
  [`der_sensitivity()`](https://joonho112.github.io/svyder/reference/der_sensitivity.md)
  and choosing a stable `tau`.
- **Selective vs blanket** correction, and why blanket adjustment is
  harmful.
- **Working with corrected draws** downstream.

> **Key insight.** The three functions form a pipeline where each stage
> annotates the same `svyder` object.
> [`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md)
> produces the raw diagnostic;
> [`der_classify()`](https://joonho112.github.io/svyder/reference/der_classify.md)
> adds tier labels and a flag column;
> [`der_correct()`](https://joonho112.github.io/svyder/reference/der_correct.md)
> adds a corrected draws matrix. Nothing is thrown away, so you can
> always step back and re-classify at a different threshold.

This vignette is scoped to the *workflow* and the *threshold*. It does
not cover how to choose the aggregation unit or variance target — see
[Choosing the Variance
Target](https://joonho112.github.io/svyder/articles/choosing-the-target.md)
for that decision — and it does not derive the tier and correction
mathematics, which live in
[`vignette("theory-classification-correction")`](https://joonho112.github.io/svyder/articles/theory-classification-correction.md).

## Step 1 — Compute

[`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md)
takes posterior draws plus the design ingredients and returns a `svyder`
object carrying the DERs and the matrices behind them. The canonical
call passes the design PSU ids as the aggregation unit:

``` r

result <- der_compute(
  nsece_demo$draws,
  y           = nsece_demo$y,
  X           = nsece_demo$X,
  group       = nsece_demo$group,
  weights     = nsece_demo$weights,
  cluster     = nsece_demo$psu,     # aggregation unit (design PSU ids)
  family      = nsece_demo$family,
  sigma_theta = nsece_demo$sigma_theta
)
```

The required inputs are:

- **`x`** — the posterior draws (an `n_draws x n_params` matrix here; S3
  methods also accept `brmsfit`, `CmdStanMCMC`, and `stanreg` fits, see
  [Backends](https://joonho112.github.io/svyder/articles/backends.md)).
- **`y`, `X`, `group`, `weights`** — the response, design matrix, group
  index, and survey weights used to fit the model.
- **`cluster`** — the aggregation unit for the sandwich meat. This is
  **required and has no default**: the DER changes materially depending
  on whether you target the design PSUs (`cluster = nsece_demo$psu`) or
  the model groups (`cluster = nsece_demo$group`). Always pass it
  explicitly.
- **`family`** — `"binomial"` here, or `"gaussian"`.
- **`sigma_theta`** — the random-effect scale (and `sigma_e` for
  gaussian models), used to form the penalized bread.

`param_types` is optional. When it is `NULL` the types are inferred from
the design matrix — a column that is constant within every group becomes
a between-group fixed effect — and a message reports the inference.
Passing `param_types` explicitly suppresses that message. On the demo
the inferred types are `intercept = fe_between`,
`poverty_cwc = fe_within`, `tiered_reim = fe_between`.

> **Note.** In the package output the poverty (within-group) covariate
> is `beta[2]` — the second design-matrix column — not `beta[1]`. The
> intercept is `beta[1]`. Keep this in mind when reading the flagged
> list below.

Printing the object *before* classification shows the diagnostic
summary:

``` r

result
#> svyder diagnostic (54 parameters)
#>   Family: binomial | N = 6785 | J = 51
#>   Target: 644 cluster(s) | meat: centered, DF-corrected | weights: unit_mean
#>   DER range: [0.235, 5.308]
#>   Tier III (DER undefined): sigma_theta
#>   (not yet classified -- run der_classify())
#>   Compute time: 0.022 sec
```

Read the header line by line. The **Target** describes the declared
variance target: 644 clusters, a meat matrix that is centered and
degrees-of-freedom corrected (the v2 defaults), and `unit_mean` weight
normalization. The **DER range** `[0.235, 5.308]` already tells you
there is real per-parameter spread — some parameters have DER well below
1 (shielded) and some well above 1 (survey-sensitive). The random-effect
scale `sigma_theta` is reported as **Tier III**: hyperparameters have no
direct data score, so their DER is undefined by construction and they
are never certified as corrected.

The returned object carries everything the later steps need:

``` r

# The parameter-level DERs
head(round(result$der, 3))
#>  beta[1]  beta[2]  beta[3] theta[1] theta[2] theta[3] 
#>    0.263    2.689    0.343    3.386    0.676    1.118

# Design metadata for the declared target
result$target[c("cluster_n", "normalize", "center", "df_correct")]
#> $cluster_n
#> [1] 644
#> 
#> $normalize
#> [1] "unit_mean"
#> 
#> $center
#> [1] TRUE
#> 
#> $df_correct
#> [1] TRUE

# Tier III parameters, excluded by construction
result$excluded
#>         param tier
#> 1 sigma_theta  III
#>                                                                                                   reason
#> 1 random-effect prior hyperparameter outside the data-level phi=(beta, theta) score block; DER undefined
```

It also stores `result$V_sand` and `result$sigma_mcmc` (the sandwich and
MCMC variances whose ratio is the DER), and `result$data` — the design
ingredients — so you can *re-target* to a different aggregation unit
without re-passing the inputs. Re-targeting is the subject of [Choosing
the Variance
Target](https://joonho112.github.io/svyder/articles/choosing-the-target.md).

## Step 2 — Classify

[`der_classify()`](https://joonho112.github.io/svyder/reference/der_classify.md)
labels each parameter with a tier and flags the ones whose DER exceeds
the threshold `tau`. The default threshold is `tau = 1.2`:

``` r

result <- der_classify(result, tau = 1.2)
#> DER Classification (tau = 1.20)
#>   Total parameters: 54
#>   Tier III excluded (DER undefined): sigma_theta
#>   Flagged: 30 (55.6%)
#>   Flagged parameters:
#>     beta[2]: DER = 2.689 [I-a] -> CORRECT
#>     theta[1]: DER = 3.386 [II] -> CORRECT
#>     theta[4]: DER = 2.210 [II] -> CORRECT
#>     theta[5]: DER = 1.568 [II] -> CORRECT
#>     theta[6]: DER = 2.104 [II] -> CORRECT
#>     theta[7]: DER = 2.232 [II] -> CORRECT
#>     theta[9]: DER = 5.308 [II] -> CORRECT
#>     theta[11]: DER = 2.655 [II] -> CORRECT
#>     theta[15]: DER = 1.576 [II] -> CORRECT
#>     theta[18]: DER = 4.025 [II] -> CORRECT
#>     theta[19]: DER = 1.992 [II] -> CORRECT
#>     theta[20]: DER = 2.469 [II] -> CORRECT
#>     theta[21]: DER = 1.790 [II] -> CORRECT
#>     theta[23]: DER = 1.311 [II] -> CORRECT
#>     theta[27]: DER = 2.324 [II] -> CORRECT
#>     theta[29]: DER = 1.631 [II] -> CORRECT
#>     theta[30]: DER = 2.968 [II] -> CORRECT
#>     theta[31]: DER = 1.288 [II] -> CORRECT
#>     theta[32]: DER = 1.715 [II] -> CORRECT
#>     theta[34]: DER = 2.844 [II] -> CORRECT
#>     theta[35]: DER = 1.258 [II] -> CORRECT
#>     theta[36]: DER = 2.280 [II] -> CORRECT
#>     theta[39]: DER = 2.222 [II] -> CORRECT
#>     theta[40]: DER = 1.567 [II] -> CORRECT
#>     theta[41]: DER = 1.821 [II] -> CORRECT
#>     theta[43]: DER = 2.301 [II] -> CORRECT
#>     theta[44]: DER = 3.217 [II] -> CORRECT
#>     theta[45]: DER = 1.505 [II] -> CORRECT
#>     theta[46]: DER = 1.225 [II] -> CORRECT
#>     theta[47]: DER = 2.134 [II] -> CORRECT
```

The four tiers partition the parameters by *how the design touches
them*:

- **Tier I-a — within-group fixed effects.** Fully exposed to the survey
  design; DER equals the Kish design effect. These are the primary
  correction candidates.
- **Tier I-b — between-group fixed effects.** Shielded by hierarchical
  shrinkage; DER equals the design effect discounted by the reliability,
  so these usually sit below the threshold.
- **Tier II — random effects.** Shrinkage-protected; DER is the design
  effect scaled by the reliability and a finite-group correction.
- **Tier III — hyperparameters.** Excluded by construction (reported in
  `$excluded`, above), never flagged.

The flag rule is a **strict** inequality: a parameter is flagged when
\mathrm{DER}\_k \> \tau, so a DER exactly equal to `tau` is *not*
flagged. On the demo this flags **30 of 54** parameters (55.6%). The
classification print shows the flagged set:

``` r

result
#> svyder diagnostic (54 parameters)
#>   Family: binomial | N = 6785 | J = 51
#>   Target: 644 cluster(s) | meat: centered, DF-corrected | weights: unit_mean
#>   DER range: [0.235, 5.308]
#>   Tier III (DER undefined): sigma_theta
#>   Threshold (tau): 1.20
#>   Flagged: 30 / 54 (55.6%)
#> 
#>   Flagged parameters:
#>     beta[2]              DER = 2.689  [I-a] -> CORRECT
#>     theta[1]             DER = 3.386  [II] -> CORRECT
#>     theta[4]             DER = 2.210  [II] -> CORRECT
#>     theta[5]             DER = 1.568  [II] -> CORRECT
#>     theta[6]             DER = 2.104  [II] -> CORRECT
#>     theta[7]             DER = 2.232  [II] -> CORRECT
#>     theta[9]             DER = 5.308  [II] -> CORRECT
#>     theta[11]            DER = 2.655  [II] -> CORRECT
#>     theta[15]            DER = 1.576  [II] -> CORRECT
#>     theta[18]            DER = 4.025  [II] -> CORRECT
#>     ... and 20 more
#>   Compute time: 0.022 sec
```

You can also work with the classification table directly. Here is the
flagged subset:

``` r

cls <- result$classification
flagged <- cls[cls$flagged, c("param_name", "param_type", "der", "tier", "action")]
head(flagged, 8)
#>    param_name param_type      der tier  action
#> 2     beta[2]  fe_within 2.689442  I-a CORRECT
#> 4    theta[1]         re 3.385886   II CORRECT
#> 7    theta[4]         re 2.209810   II CORRECT
#> 8    theta[5]         re 1.567765   II CORRECT
#> 9    theta[6]         re 2.104374   II CORRECT
#> 10   theta[7]         re 2.232233   II CORRECT
#> 12   theta[9]         re 5.307734   II CORRECT
#> 14  theta[11]         re 2.655068   II CORRECT
nrow(flagged)
#> [1] 30
```

The single flagged fixed effect is `beta[2]` (poverty, Tier I-a, DER
2.689) — a within-group covariate fully exposed to the design. The other
29 flags are random effects (`theta[j]`). The intercept `beta[1]` (DER
0.263) and the between-group policy covariate `beta[3]` (DER 0.343) are
Tier I-b and stay *below* the threshold — shrinkage shields them, so
correcting them would wrongly narrow their intervals. The tier and
correction mathematics are derived in
[`vignette("theory-classification-correction")`](https://joonho112.github.io/svyder/articles/theory-classification-correction.md).

## Step 3 — Correct

[`der_correct()`](https://joonho112.github.io/svyder/reference/der_correct.md)
rescales the flagged block toward the declared variance target and
leaves every unflagged parameter exactly as it was. The default method
is `"block_cholesky"`, the paper’s Algorithm 1:

``` r

result <- der_correct(result)
result
#> svyder diagnostic (54 parameters)
#>   Family: binomial | N = 6785 | J = 51
#>   Target: 644 cluster(s) | meat: centered, DF-corrected | weights: unit_mean
#>   DER range: [0.235, 5.308]
#>   Tier III (DER undefined): sigma_theta
#>   Threshold (tau): 1.20
#>   Flagged: 30 / 54 (55.6%)
#> 
#>   Flagged parameters:
#>     beta[2]              DER = 2.689  [I-a] -> CORRECT
#>     theta[1]             DER = 3.386  [II] -> CORRECT
#>     theta[4]             DER = 2.210  [II] -> CORRECT
#>     theta[5]             DER = 1.568  [II] -> CORRECT
#>     theta[6]             DER = 2.104  [II] -> CORRECT
#>     theta[7]             DER = 2.232  [II] -> CORRECT
#>     theta[9]             DER = 5.308  [II] -> CORRECT
#>     theta[11]            DER = 2.655  [II] -> CORRECT
#>     theta[15]            DER = 1.576  [II] -> CORRECT
#>     theta[18]            DER = 4.025  [II] -> CORRECT
#>     ... and 20 more
#> 
#>   Correction applied: 30 parameter(s) rescaled
#>   Compute time: 0.022 sec
```

The two methods differ in how they treat the *joint* structure of the
flagged block:

- **`method = "block_cholesky"`** (default, Algorithm 1) applies a joint
  Cholesky transform to the whole flagged block, matching the target
  covariance — correlations and marginals alike.
- **`method = "marginal"`** rescales each flagged parameter
  independently by \sqrt{\mathrm{DER}\_k}. It matches the target
  *marginals* but not the cross-parameter covariance. When exactly one
  parameter is flagged the two methods coincide (a scalar
  \sqrt{\mathrm{DER}} rescaling).

A third value, `method = "cholesky"`, is retired and **errors** if
requested — it was ambiguous between the block and per-parameter
transforms.

We can confirm that block-Cholesky and marginal agree when a single
parameter is flagged, using the `sim_hlr` reference dataset (which flags
exactly one parameter):

``` r

r_bc <- der_correct(
  der_classify(
    der_compute(sim_hlr$draws, y = sim_hlr$y, X = sim_hlr$X, group = sim_hlr$group,
                weights = sim_hlr$weights, cluster = sim_hlr$group,
                family = sim_hlr$family, sigma_theta = sim_hlr$sigma_theta,
                sigma_e = sim_hlr$sigma_e),
    tau = 1.2, verbose = FALSE),
  method = "block_cholesky")

r_mg <- der_correct(
  der_classify(
    der_compute(sim_hlr$draws, y = sim_hlr$y, X = sim_hlr$X, group = sim_hlr$group,
                weights = sim_hlr$weights, cluster = sim_hlr$group,
                family = sim_hlr$family, sigma_theta = sim_hlr$sigma_theta,
                sigma_e = sim_hlr$sigma_e),
    tau = 1.2, verbose = FALSE),
  method = "marginal")

sum(r_bc$classification$flagged)                       # exactly one flagged
#> [1] 1
isTRUE(all.equal(r_bc$corrected_draws, r_mg$corrected_draws))
#> [1] TRUE
```

The per-parameter scale factors are recorded in `result$scale_factors`
(an entry of `1` means the parameter was left alone):

``` r

sf <- result$scale_factors
head(round(sf, 3), 12)
#>  [1] 1.000 1.640 1.000 1.840 1.000 1.000 1.487 1.252 1.451 1.494 1.000 2.304
```

### The unflagged parameters are bitwise unchanged

This is the guarantee that makes selective correction safe: parameters
that were not flagged are copied through *bit for bit*. Every column of
the corrected matrix that corresponds to an unflagged parameter is
identical to the original draws:

``` r

flagged <- result$classification$flagged
corrected <- as.matrix(result)          # corrected draws (see below)
original  <- result$original_draws

# Unflagged columns are byte-for-byte identical
identical(corrected[, !flagged], original[, !flagged])
#> [1] TRUE

# Flagged columns did change
!identical(corrected[, flagged], original[, flagged])
#> [1] TRUE
```

The first check returns `TRUE`: nothing about the shrinkage-protected
parameters was disturbed. Only the flagged block was rescaled.

## All in one — `der_diagnose()`

[`der_diagnose()`](https://joonho112.github.io/svyder/reference/der_diagnose.md)
runs compute, classify, and correct in a single call with the same
defaults (`tau = 1.2`, `method = "block_cholesky"`, `correct = TRUE`).
It is the shortcut you would use once you have settled on a threshold:

``` r

diag_result <- der_diagnose(
  nsece_demo$draws,
  y           = nsece_demo$y,
  X           = nsece_demo$X,
  group       = nsece_demo$group,
  weights     = nsece_demo$weights,
  cluster     = nsece_demo$psu,
  family      = nsece_demo$family,
  sigma_theta = nsece_demo$sigma_theta,
  tau         = 1.2
)
```

It reaches the same classification as running the three steps by hand —
the same 30 flagged parameters:

``` r

sum(diag_result$classification$flagged)   # der_diagnose
#> [1] 30
sum(result$classification$flagged)        # step-by-step
#> [1] 30
```

Use
[`der_diagnose()`](https://joonho112.github.io/svyder/reference/der_diagnose.md)
for a quick end-to-end pass; step through
[`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md)
-\>
[`der_classify()`](https://joonho112.github.io/svyder/reference/der_classify.md)
-\>
[`der_correct()`](https://joonho112.github.io/svyder/reference/der_correct.md)
when you want to inspect the DERs before choosing `tau`, sweep several
thresholds, or correct with the `"marginal"` method.

## Threshold sensitivity

The threshold `tau` is a decision, not a law of nature, so it is worth
checking how the flagged set responds to it.
[`der_sensitivity()`](https://joonho112.github.io/svyder/reference/der_sensitivity.md)
re-classifies over a grid of thresholds (default `0.8` to `2.0`) without
recomputing the DERs:

``` r

sens <- der_sensitivity(result)
sens[, c("tau", "n_flagged", "pct_flagged")]
#>    tau n_flagged pct_flagged
#> 1  0.8        40   0.7407407
#> 2  0.9        35   0.6481481
#> 3  1.0        34   0.6296296
#> 4  1.1        34   0.6296296
#> 5  1.2        30   0.5555556
#> 6  1.3        27   0.5000000
#> 7  1.4        26   0.4814815
#> 8  1.5        26   0.4814815
#> 9  1.6        22   0.4074074
#> 10 1.7        21   0.3888889
#> 11 1.8        19   0.3518519
#> 12 1.9        18   0.3333333
#> 13 2.0        17   0.3148148
```

The flag count declines smoothly as the threshold rises:
`40, 35, 34, 34, 30` (at the default `tau = 1.2`), then
`27, 26, 26, 22, 21`, and so on. There is no cliff — the count moves by
a handful of random effects at a time — which means the classification
is not perched on a knife-edge. The default of `1.2` sits in a stable
interior region rather than at a boundary where a small change in `tau`
would flip many parameters.

> **Note.** Look for a *plateau* rather than a single “best” threshold.
> The paper reports that the correction decision is robust across \tau
> \in \[0.8, 2.0\], with 96.9% per-replication decision stability across
> the tighter band \[1.1, 1.6\] and roughly a 2% false-positive rate at
> the `1.2` default. If your flag count is stable across a band around
> your chosen `tau`, the exact value matters little.

## Selective vs blanket correction

Why correct only the flagged parameters instead of rescaling everything
by \sqrt{\mathrm{DER}}? Because a blanket adjustment applies the design
inflation to parameters that were never over-confident to begin with,
and *narrows* the ones that shrinkage was protecting.

The distinction is stark in the paper’s real-data and simulation results
(these are paper findings on the real NSECE 2019 data and its simulation
study, not outputs of the synthetic demo):

- **Real data (NSECE 2019).** A blanket correction would rescale **53 of
  54** parameters toward a rank-deficient target, with a mean
  interval-width ratio of **0.26** and a worst case of **0.04** —
  intervals collapsing to a twenty-fifth of their honest width.
  Selective correction instead widens only the parameter that needs it:
  the poverty coefficient’s 90% interval widens from **\[-0.166,
  -0.072\]** to **\[-0.207, -0.030\]**, a width ratio of **1.88**, while
  every shrinkage-protected parameter is left intact.
- **Simulation.** Under a blanket adjustment the coverage of the
  unflagged random effects collapses from about **90%** to **11–37%**.
  Selective correction keeps their coverage at the nominal level because
  it never touches them.

> **Key insight.** DER \< 1 is a warning, not a target. For a
> shrinkage-protected parameter, “correcting” toward the sandwich
> variance *narrows* the interval and destroys the hierarchical
> protection. The flagged/unflagged split exists precisely so the
> correction is applied where DER \> 1 and withheld where it would do
> harm.

The tiers and the block-Cholesky algorithm that make this split
principled are derived in
[`vignette("theory-classification-correction")`](https://joonho112.github.io/svyder/articles/theory-classification-correction.md).

## Working with corrected draws

After correction, [`as.matrix()`](https://rdrr.io/r/base/matrix.html)
returns the corrected draws if a correction has been applied, and the
original draws otherwise. This lets you write downstream code that is
agnostic to whether correction happened:

``` r

draws <- as.matrix(result)
dim(draws)
#> [1] 4000   54

# Corrected posterior summary for the flagged poverty coefficient
b2 <- draws[, "beta[2]"]
c(mean = mean(b2),
  q05  = unname(quantile(b2, 0.05)),
  q95  = unname(quantile(b2, 0.95)))
#>        mean         q05         q95 
#> -0.14944083 -0.21918324 -0.08067605
```

Feed `draws` into whatever comes next — posterior intervals,
predictions, figures, or a publication table. Because the unflagged
columns are the original draws untouched, any downstream quantity that
depends only on protected parameters is exactly what your model
produced; only the survey-sensitive parameters carry the design
correction. The
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
[`glance()`](https://generics.r-lib.org/reference/glance.html) methods
give tabular summaries of the same result and are demonstrated in the
[Case
Study](https://joonho112.github.io/svyder/articles/case-study-nsece.md).

## What’s next?

You now have the full workflow:
[`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md)
-\>
[`der_classify()`](https://joonho112.github.io/svyder/reference/der_classify.md)
-\>
[`der_correct()`](https://joonho112.github.io/svyder/reference/der_correct.md),
the
[`der_diagnose()`](https://joonho112.github.io/svyder/reference/der_diagnose.md)
shortcut, and the tools to check the threshold. Where to go from here:

- **[Choosing the Variance
  Target](https://joonho112.github.io/svyder/articles/choosing-the-target.md)**
  — the upstream decision this vignette took for granted: design-PSU vs
  model-group aggregation, weight normalization, and stratification,
  with
  [`der_compare()`](https://joonho112.github.io/svyder/reference/der_compare.md)
  to hold the two targets side by side.
- **[Case Study: NSECE-Style Survey
  Data](https://joonho112.github.io/svyder/articles/case-study-nsece.md)**
  — the whole pipeline applied end to end on the demo, with
  interpretation of the fixed and random effects, the DER decomposition,
  and a publication-ready report.
- **[`vignette("theory-classification-correction")`](https://joonho112.github.io/svyder/articles/theory-classification-correction.md)**
  — the mathematics behind the four tiers, the threshold rationale, and
  the block-Cholesky correction (Algorithm 1), including a formal
  account of why blanket correction is harmful.
