# Summarize a svyder Object

Returns a data frame with per-parameter classification details,
including tier assignment, DER value, and flagging status.

## Usage

``` r
# S3 method for class 'svyder'
summary(object, ...)
```

## Arguments

- object:

  A `svyder` object.

- ...:

  Ignored.

## Value

A `data.frame` with classification details (printed and returned
invisibly).

## See also

[`print.svyder()`](https://joonho112.github.io/svyder/reference/print.svyder.md)
for a concise summary,
[`tidy.svyder()`](https://joonho112.github.io/svyder/reference/tidy.svyder.md)
for a tidy data frame with posterior summaries.

Other svyder-methods:
[`as.matrix.svyder()`](https://joonho112.github.io/svyder/reference/as.matrix.svyder.md),
[`glance.svyder()`](https://joonho112.github.io/svyder/reference/glance.svyder.md),
[`is.svyder()`](https://joonho112.github.io/svyder/reference/is.svyder.md),
[`print.svyder()`](https://joonho112.github.io/svyder/reference/print.svyder.md),
[`tidy.svyder()`](https://joonho112.github.io/svyder/reference/tidy.svyder.md)

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
summary(result)
#>  param_name param_type       der tier                 tier_label flagged
#>     beta[1] fe_between 0.2617696  I-b        Protected (between)   FALSE
#>     beta[2]  fe_within 2.6868825  I-a           Survey-dominated    TRUE
#>     beta[3] fe_between 0.3427294  I-b        Protected (between)   FALSE
#>    theta[1]         re 3.3838218   II Protected (random effects)    TRUE
#>    theta[2]         re 0.6758803   II Protected (random effects)   FALSE
#>    theta[3]         re 1.1187639   II Protected (random effects)   FALSE
#>    theta[4]         re 2.2120536   II Protected (random effects)    TRUE
#>    theta[5]         re 1.5713900   II Protected (random effects)    TRUE
#>    theta[6]         re 2.1027473   II Protected (random effects)    TRUE
#>    theta[7]         re 2.2407594   II Protected (random effects)    TRUE
#>    theta[8]         re 0.5988151   II Protected (random effects)   FALSE
#>    theta[9]         re 5.3148375   II Protected (random effects)    TRUE
#>   theta[10]         re 0.3234469   II Protected (random effects)   FALSE
#>   theta[11]         re 2.6531533   II Protected (random effects)    TRUE
#>   theta[12]         re 0.4991804   II Protected (random effects)   FALSE
#>   theta[13]         re 0.2349819   II Protected (random effects)   FALSE
#>   theta[14]         re 0.5524033   II Protected (random effects)   FALSE
#>   theta[15]         re 1.5732575   II Protected (random effects)    TRUE
#>   theta[16]         re 0.6809155   II Protected (random effects)   FALSE
#>   theta[17]         re 0.8168582   II Protected (random effects)   FALSE
#>   theta[18]         re 4.0217038   II Protected (random effects)    TRUE
#>   theta[19]         re 1.9919767   II Protected (random effects)    TRUE
#>   theta[20]         re 2.4774293   II Protected (random effects)    TRUE
#>   theta[21]         re 1.7898408   II Protected (random effects)    TRUE
#>   theta[22]         re 0.8740200   II Protected (random effects)   FALSE
#>   theta[23]         re 1.3111319   II Protected (random effects)    TRUE
#>   theta[24]         re 0.9736441   II Protected (random effects)   FALSE
#>   theta[25]         re 0.7492859   II Protected (random effects)   FALSE
#>   theta[26]         re 1.1561149   II Protected (random effects)   FALSE
#>   theta[27]         re 2.3264034   II Protected (random effects)    TRUE
#>   theta[28]         re 1.1845663   II Protected (random effects)   FALSE
#>   theta[29]         re 1.6319237   II Protected (random effects)    TRUE
#>   theta[30]         re 2.9722557   II Protected (random effects)    TRUE
#>   theta[31]         re 1.2897827   II Protected (random effects)    TRUE
#>   theta[32]         re 1.7132695   II Protected (random effects)    TRUE
#>   theta[33]         re 0.7049967   II Protected (random effects)   FALSE
#>   theta[34]         re 2.8424594   II Protected (random effects)    TRUE
#>   theta[35]         re 1.2585246   II Protected (random effects)    TRUE
#>   theta[36]         re 2.2770740   II Protected (random effects)    TRUE
#>   theta[37]         re 0.7444208   II Protected (random effects)   FALSE
#>   theta[38]         re 1.1215973   II Protected (random effects)   FALSE
#>   theta[39]         re 2.2220272   II Protected (random effects)    TRUE
#>   theta[40]         re 1.5670170   II Protected (random effects)    TRUE
#>   theta[41]         re 1.8183090   II Protected (random effects)    TRUE
#>   theta[42]         re 0.8853889   II Protected (random effects)   FALSE
#>   theta[43]         re 2.2998027   II Protected (random effects)    TRUE
#>   theta[44]         re 3.2334813   II Protected (random effects)    TRUE
#>   theta[45]         re 1.5147483   II Protected (random effects)    TRUE
#>   theta[46]         re 1.2276126   II Protected (random effects)    TRUE
#>   theta[47]         re 2.1405640   II Protected (random effects)    TRUE
#>   theta[48]         re 0.8888311   II Protected (random effects)   FALSE
#>   theta[49]         re 0.7703205   II Protected (random effects)   FALSE
#>   theta[50]         re 0.8947901   II Protected (random effects)   FALSE
#>   theta[51]         re 0.7644108   II Protected (random effects)   FALSE
#>   action
#>   retain
#>  CORRECT
#>   retain
#>  CORRECT
#>   retain
#>   retain
#>  CORRECT
#>  CORRECT
#>  CORRECT
#>  CORRECT
#>   retain
#>  CORRECT
#>   retain
#>  CORRECT
#>   retain
#>   retain
#>   retain
#>  CORRECT
#>   retain
#>   retain
#>  CORRECT
#>  CORRECT
#>  CORRECT
#>  CORRECT
#>   retain
#>  CORRECT
#>   retain
#>   retain
#>   retain
#>  CORRECT
#>   retain
#>  CORRECT
#>  CORRECT
#>  CORRECT
#>  CORRECT
#>   retain
#>  CORRECT
#>  CORRECT
#>  CORRECT
#>   retain
#>   retain
#>  CORRECT
#>  CORRECT
#>  CORRECT
#>   retain
#>  CORRECT
#>  CORRECT
#>  CORRECT
#>  CORRECT
#>  CORRECT
#>   retain
#>   retain
#>   retain
#>   retain
```
