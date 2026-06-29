# Print summary of mixture fit

Displays formatted summary statistics for a fitted mixture model.

## Usage

``` r
# S3 method for class 'summary.shrinkr_mixture'
print(x, digits = 3, ...)
```

## Arguments

- x:

  A summary object from
  [`summary.shrinkr_mixture`](summary.shrinkr_mixture.md).

- digits:

  Number of decimal digits to display. Default is 3.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the input object `x`.

## See also

[`summary.shrinkr_mixture`](summary.shrinkr_mixture.md) for creating
summaries

## Examples

``` r
set.seed(1)
samples <- list(
  group1 = matrix(rnorm(100, 0.0, 0.5), ncol = 1),
  group2 = matrix(rnorm(100, 0.5, 0.5), ncol = 1)
)
mix <- fit_mixture(samples, K_max = 2, verbose = FALSE)
summ <- summary(mix)
print(summ)
#> == Mixture Model Summary ============================
#> 
#> Model: XII 
#> Components: 1 
#> Log-likelihood: -129.32 
#> BIC: -272.456 
#> Parameters: 3 
#> 
#> Component weights:
#> # A tibble: 1 × 2
#>   component weight
#>       <int>  <dbl>
#> 1         1      1
#> 
#> Variable-wise summaries (weighted across components):
#> # A tibble: 2 × 5
#>   variable weighted_mean weighted_sd range_mean range_sd
#>   <chr>            <dbl>       <dbl>      <dbl>    <dbl>
#> 1 group1          0.0544       0.462          0        0
#> 2 group2          0.481        0.462          0        0
#> 
#> -----------------------------------------------------
```
