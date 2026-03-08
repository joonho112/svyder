# Create a ggplot2 Visualization of DER Results

Provides the same functionality as
[`plot.svyder()`](https://joonho112.github.io/svyder/reference/plot.svyder.md)
but requires ggplot2 and always returns a `ggplot` object. This is the
preferred method when using ggplot2 directly.

## Usage

``` r
autoplot.svyder(
  object,
  type = c("profile", "decomposition", "comparison"),
  ...
)
```

## Arguments

- object:

  A `svyder` object.

- type:

  Character; plot type: `"profile"` (default), `"decomposition"`, or
  `"comparison"`.

- ...:

  Additional arguments passed to the underlying plot function.

## Value

A `ggplot` object.

## See also

[`plot.svyder()`](https://joonho112.github.io/svyder/reference/plot.svyder.md)
for the generic plot method.

Other visualization:
[`plot.svyder()`](https://joonho112.github.io/svyder/reference/plot.svyder.md)

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
  ggplot2::autoplot(result, type = "profile")
}

```
