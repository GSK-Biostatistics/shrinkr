# Summary statistics for mixture fits

Computes comprehensive summary statistics for fitted mixture components,
including component-wise and variable-wise summaries.

## Usage

``` r
# S3 method for class 'shrinkr_mixture'
summary(object, ...)
```

## Arguments

- object:

  A `shrinkr_mixture` object from [`fit_mixture`](fit_mixture.md).

- ...:

  Additional arguments (currently unused).

## Value

A list with class "summary.shrinkr_mixture" containing:

- components:

  Data frame with component weights and sizes

- by_variable:

  Data frame with per-variable mixture summaries

- model_info:

  List with model selection details

- diagnostics:

  Fit diagnostics

## See also

[`fit_mixture`](fit_mixture.md) for fitting mixture models,
[`print.shrinkr_mixture`](print.shrinkr_mixture.md) for quick overview

## Examples

``` r
set.seed(1)
samples <- list(
  group1 = matrix(rnorm(100, 0.0, 0.5), ncol = 1),
  group2 = matrix(rnorm(100, 0.5, 0.5), ncol = 1)
)
mix <- fit_mixture(samples, K_max = 2, verbose = FALSE)
summ <- summary(mix)
summ$components
#> # A tibble: 1 × 2
#>   component weight
#>       <int>  <dbl>
#> 1         1      1
```
