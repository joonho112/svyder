# Plot DER Diagnostic Results

Creates diagnostic visualizations for Design Effect Ratio results. Three
plot types are available:

- profile:

  DER values by parameter with tier-based coloring and threshold line.
  Default plot type.

- decomposition:

  Scatter of DER against decomposition components (shrinkage factor B,
  design effect DEFF).

- comparison:

  Side-by-side naive vs corrected credible intervals. Requires
  [`der_correct()`](https://joonho112.github.io/svyder/reference/der_correct.md)
  to have been applied.

## Usage

``` r
# S3 method for class 'svyder'
plot(x, type = c("profile", "decomposition", "comparison"), ...)
```

## Arguments

- x:

  A `svyder` object.

- type:

  Character; plot type: `"profile"` (default), `"decomposition"`, or
  `"comparison"`.

- ...:

  Additional arguments passed to the underlying plot function.

## Value

If ggplot2 is available, a `ggplot` object (invisibly). Otherwise, `x`
is returned invisibly.

## Details

If ggplot2 is installed, returns a ggplot object. Otherwise, falls back
to base R graphics.

## See also

[`autoplot.svyder()`](https://joonho112.github.io/svyder/reference/autoplot.svyder.md)
for explicit ggplot2 usage.

Other visualization:
[`autoplot.svyder()`](https://joonho112.github.io/svyder/reference/autoplot.svyder.md)

## Examples

``` r
data(nsece_demo)
result <- der_diagnose(
  nsece_demo$draws,
  y = nsece_demo$y, X = nsece_demo$X,
  group = nsece_demo$group, weights = nsece_demo$weights,
  psu = nsece_demo$psu, family = "binomial",
  sigma_theta = nsece_demo$sigma_theta,
  param_types = nsece_demo$param_types
)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot(result, type = "profile")
}

```
