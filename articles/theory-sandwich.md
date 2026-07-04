# The Declared Variance Target: Bread and Meat

## Overview

Every Design Effect Ratio is a ratio, and its numerator is a single
object: the **declared variance target** V\_{\mathrm{target}} = H^{-1}
J_c H^{-1}, a design-based sandwich covariance evaluated at the
posterior mean \hat{\boldsymbol{\phi}}. This article derives each factor
— the *bread* H and the *meat* J_c — and, more importantly, explains the
modelling choices that promote V\_{\mathrm{target}} from an ad-hoc
formula to a **well-defined estimand**: which Hessian, which aggregation
unit, which weight scaling, and how the cluster universe is declared.

This vignette is written for methodologists who want to know precisely
what svyder computes before it takes a ratio. We separate what is
established from what is new:

- **Established.** The design-based sandwich form H^{-1} J_c H^{-1} is
  the survey-statistics standard: linearization variance in the style of
  Binder (1983), with the finite-cluster degrees-of-freedom correction
  that underlies survey-package variance estimators. The
  pseudo-posterior it corrects is the weighted-likelihood construction
  of Savitsky & Williams (2022), and the posterior covariance
  calibration is the Williams & Savitsky (2021) covariance-calibration
  transform.
- **Novel (v0.2.0).** The *declared-target* formalization — treating the
  aggregation unit, weight convention, and cluster universe as explicit,
  recorded parts of the estimand rather than silent defaults — together
  with the `cluster_strata` handling of selected-but-empty second-stage
  clusters and a corrected gaussian score.

> **Key insight.** V\_{\mathrm{target}} is not “the” design variance. It
> is *a* design variance, indexed by choices you must declare. Changing
> the aggregation unit or the weight scaling changes the estimand, not
> merely the estimate. svyder records these choices in the object’s
> `$target` slot so that a DER is always interpretable.

We cover, in order: the bread H and why it must be the *posterior*
Hessian; the meat J_c with stratum centering and the Binder correction;
the aggregation unit as part of the estimand; weight conventions and
scale sensitivity; the cluster universe and empty PSUs; and the v0.2.0
gaussian score fix. The downstream question — what these choices *imply*
for DER values in closed form — is the subject of
[`vignette("theory-decomposition")`](https://joonho112.github.io/svyder/articles/theory-decomposition.md);
here we build the target, we do not decompose it.

We use the paper’s notation throughout. The parameter vector is
\boldsymbol{\phi} = (\boldsymbol{\beta}^{\mathsf{T}},
\boldsymbol{\theta}^{\mathsf{T}})^{\mathsf{T}} of dimension d = p + J,
stacking p fixed effects and J group-level random intercepts;
\tilde{w}\_{ij} is the normalized weight for unit i in group j; and
f(y\_{ij} \mid \mathbf{x}\_{ij}^{\mathsf{T}}\boldsymbol{\beta} +
\theta_j) is the observation-level likelihood.

## The bread H: penalized posterior Hessian

The bread is the negative curvature of the log-*objective* that the
sampler explores — the penalized pseudo-log-likelihood — evaluated at
the posterior mean: H \\=\\ -\sum\_{j,i} \tilde{w}\_{ij}\\
\nabla^2\_{\boldsymbol{\phi}} \log f\\\left(y\_{ij} \mid
\mathbf{x}\_{ij}^{\mathsf{T}}\boldsymbol{\beta} + \theta_j\right)
\bigg\|\_{\hat{\boldsymbol{\phi}}} \\+\\
\operatorname{diag}\\\left(\mathbf{0}\_p, \tau, \ldots, \tau\right),
\qquad \tau = \frac{1}{\sigma\_\theta^2}.

Two terms build it. The first is the **weighted likelihood curvature**,
summed over all sampled units and contributed by the data layer alone.
The second is the **random-effect prior precision** \tau, one entry on
each of the J diagonal positions of the \boldsymbol{\theta} block. The
fixed-effect prior contributes nothing, by the
diffuse-\boldsymbol{\beta} convention discussed below.

### Why the *posterior* Hessian, not the likelihood Hessian

The choice of H is not cosmetic. The pure pseudo-likelihood Hessian —
the first term on its own — is **singular** for any hierarchical model
that carries both an intercept and a full set of group effects. The
reason is a flat direction. Consider perturbing the intercept up by a
constant while shifting every group effect down by the same constant:
(\beta_0, \boldsymbol{\theta}) \\\longmapsto\\ (\beta_0 + \varepsilon,\\
\boldsymbol{\theta} - \varepsilon \mathbf{1}\_J). The linear predictor
\mathbf{x}\_{ij}^{\mathsf{T}}\boldsymbol{\beta} + \theta_j is unchanged
for every unit, so the log-likelihood is exactly flat along the
direction (1, \mathbf{0}\_{p-1}, -\mathbf{1}\_J). The likelihood Hessian
therefore has this vector in its null space and cannot be inverted — yet
the sandwich needs H^{-1}.

The random-effect prior precision removes *exactly* this flat direction.
The term \operatorname{diag}(\mathbf{0}\_p, \tau, \ldots, \tau) places
curvature \tau \> 0 on the \boldsymbol{\theta} coordinates, and the flat
direction has nonzero \boldsymbol{\theta}-component, so it is no longer
in the null space of H. The penalized Hessian is positive definite, and
H^{-1} exists.

> **Note.** This is why svyder inverts the *posterior* Hessian and not
> the likelihood Hessian. The bread must be the curvature of the
> objective the MCMC sampler actually explores. Using the likelihood
> Hessian would not only be a different estimand — it would not be
> invertible.

### The diffuse-\boldsymbol{\beta} convention

By default the fixed effects carry a flat (improper-limit) prior, so
they add **no curvature** to the bread. In the API this is
`beta_prior_sd = Inf`, the diffuse-\boldsymbol{\beta} convention. It
matches the modelling assumption under which the closed-form
decompositions hold (flat prior on \boldsymbol{\beta}; see
[`vignette("theory-decomposition")`](https://joonho112.github.io/svyder/articles/theory-decomposition.md))
and keeps the fixed-effect block of H equal to the weighted Fisher
information alone.

Supplying a finite value instead reproduces an
**informative-\boldsymbol{\beta} bread**: a term 1/\sigma\_\beta^2
enters the fixed-effect diagonal, shrinking the fixed effects toward
zero and shifting the target accordingly. The default of \infty is the
neutral choice; a finite value is a deliberate declaration that the
prior on \boldsymbol{\beta} is part of the target.

### \tau is a post-fit plug-in, not a circularity

The precision \tau = 1/\sigma\_\theta^2 depends on the random-effect
variance, which is itself estimated. svyder evaluates it at the
posterior plug-in \hat{\tau} = 1/\hat{\sigma}\_\theta^2. This does
**not** make the diagnostic circular. The DER is a *conditional,
post-fit* diagnostic: it asks, given the fitted posterior — including
the fitted \hat{\sigma}\_\theta^2 — how the model’s reported uncertainty
compares with the design-based target. Uncertainty *in*
\hat{\sigma}\_\theta^2 is a separate question, and hyperparameters are
Tier III. For \sigma\_\theta, the score belongs to the unweighted
random-effect prior rather than the observation-level \boldsymbol\phi
score. For gaussian \sigma_e, svyder conditions on the plug-in
dispersion and constructs only the location-block sandwich. We return to
that exclusion when the excluded set appears in the computed object
below; the tier scheme itself is the subject of
[`vignette("theory-classification-correction")`](https://joonho112.github.io/svyder/articles/theory-classification-correction.md).

In the demo object, the plug-in precision comes from
`sigma_theta = 0.66`, and svyder records the value it used:

``` r

obj <- der_compute(
  nsece_demo$draws,
  y           = nsece_demo$y,
  X           = nsece_demo$X,
  group       = nsece_demo$group,
  weights     = nsece_demo$weights,
  cluster     = nsece_demo$psu,        # design PSUs -- the aggregation unit
  family      = nsece_demo$family,
  sigma_theta = nsece_demo$sigma_theta,
  param_types = nsece_demo$param_types
)

# The bread is a full d x d matrix (d = p + J = 3 + 51 = 54 here)
dim(obj$H_obs)
#> [1] 54 54
# svyder records the plug-in convention and the sigma_theta it used
obj$target$plug_in
#> [1] "posterior_mean"
obj$target$beta_prior_sd
#> [1] Inf
obj$target$sigma_theta
#> [1] 0.66
```

## The meat J_c: cluster-aggregated score

The meat is the design-based variability of the *estimating function*.
It is built by aggregating the weighted score to the cluster level,
centering within strata, and applying a finite-cluster correction: J_c
\\=\\ \sum_h \frac{C_h}{C_h - 1} \sum\_{c \in h}
\left(\mathbf{t}\_{hc} - \bar{\mathbf{t}}\_h\right)
\left(\mathbf{t}\_{hc} - \bar{\mathbf{t}}\_h\right)^{\mathsf{T}}, where
the **cluster score total** for cluster c in stratum h collects the
weighted gradient over every unit the cluster contains: \mathbf{t}\_{hc}
\\=\\ \sum\_{(i,j)\in c} \tilde{w}\_{ij}\\ \nabla\_{\boldsymbol{\phi}}
\log f\\\left(y\_{ij} \mid
\mathbf{x}\_{ij}^{\mathsf{T}}\boldsymbol{\beta} + \theta_j\right)
\bigg\|\_{\hat{\boldsymbol{\phi}}}.

This is the survey-linearization meat: a between-cluster sum of squares
of the score totals, with clusters as the ultimate sampling units.

### Stratum centering

The centered contribution uses \mathbf{t}\_{hc} - \bar{\mathbf{t}}\_h,
where \bar{\mathbf{t}}\_h = C_h^{-1} \sum\_{c \in h} \mathbf{t}\_{hc} is
the stratum mean of the cluster score totals. Centering is not optional
bookkeeping. The weighted score totals **do not vanish at the posterior
mean**: the posterior mean is a Bayesian summary, not a weighted
maximum-likelihood root, so the sample estimating equation is not solved
exactly at \hat{\boldsymbol{\phi}}. A design-based variance of an
estimating function must be centered at its (stratum-wise) mean, exactly
as a survey variance estimator subtracts the stratum mean of the PSU
totals. In the API this is `center_meat = TRUE`, and setting it to
`FALSE` reverts to the uncentered form \sum_c \mathbf{t}\_c
\mathbf{t}\_c^{\mathsf{T}} that assumes the score sum is zero —
appropriate only if you are matching a maximum-likelihood convention
exactly.

### The Binder finite-cluster correction

Each stratum carries the factor C_h / (C_h - 1), the standard Binder
(1983) degrees-of-freedom adjustment for a stratum containing C_h
clusters. It corrects the downward bias of a between-cluster sum of
squares built from a finite number of clusters — the same C_h/(C_h-1)
that survey-package variance estimators apply. It matters most where it
bites hardest: small-cluster strata, in which the naive sum of squares
would materially understate the variance. In the API this is
`df_correct = TRUE`. Together,
`center_meat = FALSE, df_correct = FALSE, strata = NULL` reproduce the
uncorrected v1 meat.

> **Note.** The two flags carry distinct statistical content.
> `center_meat` decides *what deviations* are squared (from the stratum
> mean, or from zero); `df_correct` decides *how they are scaled* (by
> C_h/(C_h-1), or not). svyder records both in `$target` so the meat’s
> construction is never ambiguous.

``` r

# The recorded meat convention for the computed object
obj$target$center       # stratum centering on?
#> [1] TRUE
obj$target$df_correct   # Binder C_h/(C_h - 1) on?
#> [1] TRUE
obj$target$strata_n     # number of strata (1 = single-stratum meat here)
#> [1] 1
```

### Rank, and why the diagonal survives

The meat is a sum of outer products, one (centered) term per cluster, so
\operatorname{rank}(J_c) \\\le\\ \\\\\text{clusters}\\. When there are
fewer clusters than parameters d, or when centering removes a degree of
freedom, J_c is **rank-deficient** and the full sandwich
V\_{\mathrm{target}} is singular. That is not a problem for the DER. The
DER reads only the **diagonal** entries \[V\_{\mathrm{target}}\]\_{kk},
and those stay well-defined even when the matrix as a whole is not
positive definite. Rank deficiency does constrain *joint* corrections
spanning many flagged parameters — svyder falls back to a nearest
positive-definite target there — but that is a correction-step concern,
taken up in
[`vignette("theory-classification-correction")`](https://joonho112.github.io/svyder/articles/theory-classification-correction.md).

``` r

# One (centered) outer product per cluster caps the meat's rank
dim(obj$J_cluster)        # d x d = 54 x 54
#> [1] 54 54
qr(obj$J_cluster)$rank    # bounded by min(#clusters, d); diagonal still usable
#> [1] 52
```

## The aggregation unit is part of the estimand

The cluster index c in J_c is a *modelling choice*, and it is the single
most consequential one. The formula does not change — the same centering
and the same C_h/(C_h-1) apply — but *what counts as a cluster*
determines which question the target answers.

- **Design-PSU target** (`cluster = <design PSU ids>`). Clusters are the
  survey’s own primary sampling units. The meat then captures the
  design’s actual sampling variability: the variance an analyst would
  attribute to the sampling mechanism that generated the data.
- **Model-group target** (`cluster = group`). Clusters are the model’s
  random-effect groups. The meat captures variability of the group-level
  structure the model posits.

These are different estimands and they frequently disagree. On the demo,
the design-PSU target flags 30 of 54 parameters at \tau = 1.2; the
model-group target flags 25 of 54. Same data, same model, different
aggregation unit, different diagnostic answer. This is why svyder makes
`cluster` a **required** argument with no default: the package refuses
to guess an estimand on your behalf.

> **Key insight.** There is no “correct” aggregation unit in the
> abstract. The design-PSU target answers “does the posterior reflect
> the sampling design?”; the model-group target answers “does the
> posterior reflect its own group structure?”. Choosing between them is
> a substantive decision, treated in the applied companion [Choosing the
> Variance
> Target](https://joonho112.github.io/svyder/articles/choosing-the-target.md).

svyder records which unit was used, so the estimand is legible from the
object:

``` r

obj$target$cluster_is_group   # FALSE -> design-PSU target (we passed cluster = psu)
#> [1] FALSE
obj$target$cluster_n          # 644 design PSUs in the demo
#> [1] 644
```

## Weight conventions

The normalized weight \tilde{w}\_{ij} enters *both* the pseudo-posterior
and the target, so **how the weights are scaled is itself part of the
estimand.** svyder supports three conventions via `normalize`, and
records the choice in `$target$normalize`.

- **Global unit-mean** (`"unit_mean"`, the default and the paper’s
  primary convention): \tilde{w}\_{ij} \\=\\ N \cdot
  \frac{w\_{ij}}{\sum\_{i',j'} w\_{i',j'}}, which scales the whole
  sample to mean weight 1 and **preserves between-group weight
  variation**.
- **Within-group** (`"group_size"`): \sum_i \tilde{w}\_{ij} = n_j \quad
  \text{for every group } j, which forces each group’s weights to sum to
  its own sample size and thereby **removes between-group weight
  variation**, reallocating effective sample sizes across groups.
- **None** (`"none"`): the raw weights are passed through unscaled.

For fixed effects the convention is often close to immaterial, but
random-effect DERs are **not scale-invariant**: they depend on
group-level effective sample sizes, which the two conventions apportion
differently. The paper’s NSECE analysis shows the convention can
*reverse* the ordering of the two targets — the poverty coefficient’s
model-group DER is 4.59 under unit-mean but 2.70 under within-group,
while its design-PSU DER moves the other way, from 3.55 to 4.70. Because
it changes the answer, the weight convention must be declared alongside
the aggregation unit.

> **Note.** The bundled `nsece_demo` is deliberately *not* a good place
> to see this: its weights are already within-group normalized with mean
> 1, so `"unit_mean"`, `"group_size"`, and `"none"` coincide. The
> weight-convention sensitivity is a fact about the real NSECE data
> (paper, Section 5.4), not something the demo can exhibit. The applied
> [Choosing the Variance
> Target](https://joonho112.github.io/svyder/articles/choosing-the-target.md)
> vignette discusses the practical choice.

``` r

obj$target$normalize          # "unit_mean" -- the declared scaling
#> [1] "unit_mean"
obj$target$weights_raw_mean   # mean of the raw weights before normalization
#> [1] 1
```

## The cluster universe and empty PSUs

A two-stage design selects first-stage clusters and then, within each, a
second stage. A subtlety arises when a selected cluster contributes **no
observations** to the analysis sample — for instance a PSU that was in
the design but whose sampled units all fall outside the modelled
subpopulation. If such a cluster simply disappears from the data, it
silently drops out of both the centering and the degrees-of-freedom
count, and the meat is computed as though the design had fewer clusters
than it truly did.

The `cluster_strata` argument fixes this by **declaring the full cluster
universe**. It supplies the stratum label of every cluster that *could*
have contributed, whether or not it did. Selected-but-empty clusters
then enter J_c as **zero score-total clusters**: their \mathbf{t}\_{hc}
= \mathbf{0}, so they add no outer product of their own, but they still

- shift the stratum mean \bar{\mathbf{t}}\_h used for centering (a zero
  total is a real, informative total, not a missing one), and
- raise C_h, so they participate in the C_h/(C_h-1) correction.

Both effects make the meat reflect the design that was actually run
rather than only the clusters that happened to yield data.

> **Note.** `cluster_strata` requires integer cluster ids running over
> `1:length(cluster_strata)` and requires `strata` to be supplied,
> because an empty cluster has no data from which to infer its stratum.
> Leave it `NULL` — as on the demo, where every realized cluster has at
> least one observation — and no zero-total clusters are added. The demo
> object confirms the empty count is zero and the universe was not
> separately declared:

``` r

obj$target$cluster_n_empty            # 0 empty clusters in the demo
#> [1] 0
obj$target$cluster_universe_declared  # FALSE -- no cluster_strata supplied here
#> [1] FALSE
```

This zero-total handling is the second piece of the v0.2.0
declared-target machinery: it lets you specify the *population* of
clusters, not merely the realized ones, so the design-based meat is
faithful to the sampling plan.

## The gaussian score fix (v0.2.0)

The final v0.2.0 correction concerns the *gaussian* family. The score
that enters both \mathbf{t}\_{hc} (the meat) and the curvature (the
bread) must be scaled consistently, and for gaussian likelihoods v1 was
not.

The correct weighted gaussian score residual is \nabla\_\phi \log
f(y\_{ij} \mid \mu\_{ij})\big\|\_{\hat{\phi}} \\=\\
\tilde{w}\_{ij}\\\frac{y\_{ij} - \mu\_{ij}}{\sigma_e^2}, so that the
cluster score totals become \mathbf{t}\_{hc} \\=\\ \sum\_{(i,j)\in c}
\tilde{w}\_{ij}\\ \frac{y\_{ij} - \mu\_{ij}}{\sigma_e^2}. The
1/\sigma_e^2 factor is the gaussian **working weight**: in IRLS terms
the effective weight of a gaussian observation is
\tilde{w}\_{ij}/\sigma_e^2, and the score must carry the same
\sigma_e^{-2} that the working weight does so that bread and meat are
scaled on the same footing.

Version 1 used \tilde{w}\_{ij}(y\_{ij} - \mu\_{ij}) — dropping the
\sigma_e^{-2}. Because the score entered the meat quadratically, this
made the gaussian DER scale like \sigma_e^{4}, which is wrong unless
\sigma_e = 1. The error was invisible whenever the error scale happened
to be one and grew with \sigma_e otherwise. The **binomial family is
unaffected**: its natural score \tilde{w}\_{ij}(y\_{ij} - \mu\_{ij})
carries no separate dispersion parameter, so the v1 form was already
correct for it. Non-gaussian families likewise use their natural score;
only the gaussian case needed the working-weight correction, and with it
the gaussian target now aligns exactly with survey-weighted linear
regression variance (Binder, 1983).

To ground the discussion, the `$target` slot of a computed object
exposes the full metadata of the estimand — every declaration described
in this vignette, in one place:

``` r

# The complete declared-target record for the demo object
str(obj$target)
#> List of 12
#>  $ cluster_n                : int 644
#>  $ cluster_n_empty          : int 0
#>  $ cluster_universe_declared: logi FALSE
#>  $ cluster_is_group         : logi FALSE
#>  $ strata_n                 : int 1
#>  $ center                   : logi TRUE
#>  $ df_correct               : logi TRUE
#>  $ normalize                : chr "unit_mean"
#>  $ weights_raw_mean         : num 1
#>  $ plug_in                  : chr "posterior_mean"
#>  $ beta_prior_sd            : num Inf
#>  $ sigma_theta              : num 0.66
```

For a gaussian fit the same machinery applies, and the corrected score
is what
[`der_compute()`](https://joonho112.github.io/svyder/reference/der_compute.md)
uses under the hood. The bundled `sim_hlr` object is gaussian with
\sigma_e = 1, so it is exactly the case in which the v1 and v0.2.0
scores *coincide* — a useful reference precisely because the fix leaves
it unchanged. Note that a gaussian model has a second variance
component, \sigma_e. svyder uses it in the gaussian working weights
above, but still excludes \sigma_e as Tier III because the declared DER
target is conditional on that plug-in dispersion and covers only the
location block \boldsymbol\phi:

``` r

g <- der_compute(
  sim_hlr$draws,
  y           = sim_hlr$y,
  X           = sim_hlr$X,
  group       = sim_hlr$group,
  weights     = sim_hlr$weights,
  cluster     = sim_hlr$group,       # model-group target for this reference fit
  family      = sim_hlr$family,      # "gaussian"
  sigma_theta = sim_hlr$sigma_theta,
  sigma_e     = sim_hlr$sigma_e,     # 1.0 here -> v1 and v2 scores agree
  param_types = sim_hlr$param_types
)

# Both variance components are excluded from the DER (Tier III)
g$excluded
#>         param tier
#> 1 sigma_theta  III
#> 2     sigma_e  III
#>                                                                                                                                      reason
#> 1                                    random-effect prior hyperparameter outside the data-level phi=(beta, theta) score block; DER undefined
#> 2 gaussian dispersion plug-in; svyder conditions on sigma_e and constructs the location-block sandwich for phi=(beta, theta); DER undefined
```

> **Key insight.** The gaussian fix is a statement about *consistency*,
> not magnitude. Bread and meat must be built from the same working
> weights, or the sandwich no longer targets the survey-weighted
> regression variance it claims to. With \sigma_e = 1 nothing changes;
> the fix earns its keep the moment the error scale departs from one.

## What’s next?

You now have the full construction of the numerator that every DER
divides by: the penalized-posterior bread H, the cluster-aggregated meat
J_c with stratum centering and the Binder correction, and the
aggregation, weight, and cluster-universe declarations that make
V\_{\mathrm{target}} = H^{-1} J_c H^{-1} a well-defined estimand rather
than a formula with hidden defaults.

- To see what these choices *imply* — the closed-form expressions that
  turn this target into DER values, the fixed- and random-effect
  decompositions, and the conservation law — read
  [`vignette("theory-decomposition")`](https://joonho112.github.io/svyder/articles/theory-decomposition.md).
- To make the aggregation-unit and weight-convention choices in
  practice, and to compare targets on real numbers, see the applied
  [Choosing the Variance
  Target](https://joonho112.github.io/svyder/articles/choosing-the-target.md).

## References

Binder, D. A. (1983). On the variances of asymptotically normal
estimators from complex surveys. *International Statistical Review*,
51(3), 279–292.

Kish, L. (1965). *Survey Sampling*. Wiley.

Lee, J., Williams, M. R., & Savitsky, T. D. (2026). Design Effect Ratios
for Bayesian Survey Models: A Diagnostic Framework for Identifying
Survey-Sensitive Parameters. *Journal of Survey Statistics and
Methodology* (submitted).

Savitsky, T. D., & Williams, M. R. (2022). Pseudo Bayesian mixed models
under informative sampling. *Journal of Official Statistics*, *38*(4),
901–928. <https://doi.org/10.2478/jos-2022-0039>

Williams, M. R., & Savitsky, T. D. (2021). Uncertainty estimation for
pseudo-Bayesian inference under complex sampling. *International
Statistical Review*, *89*(1), 72–107.
<https://doi.org/10.1111/insr.12376>
