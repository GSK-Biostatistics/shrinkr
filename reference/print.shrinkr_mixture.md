# Print method for mixture fits

Displays summary information about a fitted mixture model in a readable
format.

## Usage

``` r
# S3 method for class 'shrinkr_mixture'
print(x, ...)
```

## Arguments

- x:

  A `shrinkr_mixture` object from [`fit_mixture`](fit_mixture.md).

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the input object `x`.

## See also

[`fit_mixture`](fit_mixture.md) for fitting mixture models,
[`summary.shrinkr_mixture`](summary.shrinkr_mixture.md) for detailed
summaries

## Examples

``` r
set.seed(1)
samples <- list(
  group1 = matrix(rnorm(100, 0.0, 0.5), ncol = 1),
  group2 = matrix(rnorm(100, 0.5, 0.5), ncol = 1)
)
mix <- fit_mixture(samples, K_max = 2, verbose = FALSE)
print(mix)
#> ── Gaussian Mixture Model Fit ──────────────────────────────────────────────────
#> 
#> Components: 1
#> Variables: 2 (group1, group2)
#> Samples: 100
#> 
#> Model: XII
#> BIC: -272.46
#> 
#> 
#> ── Component weights 
#> k = 1: 1
#> 
#> ────────────────────────────────────────────────────────────────────────────────
#> ℹ Use `plot()` to visualize marginal fits
#> ℹ Use `as.data.frame()` for component specifications
#> ℹ Use `summary()` for detailed statistics
```
