# Print summary of prior predictive samples

Displays formatted summary statistics for prior predictive samples.

## Usage

``` r
# S3 method for class 'summary.shrinkr_prior_pred'
print(x, digits = 3, ...)
```

## Arguments

- x:

  A summary object from
  [`summary.shrinkr_prior_pred`](summary.shrinkr_prior_pred.md).

- digits:

  Number of decimal digits to display. Default is 3.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the input object `x`.

## See also

[`summary.shrinkr_prior_pred`](summary.shrinkr_prior_pred.md) for
creating summaries

## Examples

``` r
priors <- list(
  mu = distributional::dist_normal(0, 5),
  tau = distributional::dist_truncated(distributional::dist_normal(0, 1), lower = 0)
)
prior_pred <- sample_prior_predictive(priors, n_groups = 3, n_draws = 50)
print(summary(prior_pred))
#> == Prior Predictive Summary =========================
#> 
#> Based on 50 draws for 3 groups
#> 
#> Hyperparameters:
#> # A tibble: 2 × 6
#>   parameter   mean    sd     q2.5  q50.0 q97.5
#>   <chr>      <dbl> <dbl>    <dbl>  <dbl> <dbl>
#> 1 mu        -0.157 5.47  -11.1    -0.145  8.89
#> 2 tau        0.879 0.658   0.0407  0.769  2.61
#> 
#> Subgroup effects (theta):
#> # A tibble: 3 × 6
#>   group     mean    sd  q2.5   q50.0 q97.5
#>   <chr>    <dbl> <dbl> <dbl>   <dbl> <dbl>
#> 1 group1 -0.209   5.57 -10.7 -0.341  10.8 
#> 2 group2 -0.0604  5.64 -10.9  0.190  10.0 
#> 3 group3 -0.0191  5.70 -10.7 -0.0812  9.51
#> 
#> -----------------------------------------------------
```
