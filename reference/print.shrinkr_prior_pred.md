# Print method for prior predictive samples

Displays summary information about prior predictive samples in a
readable format.

## Usage

``` r
# S3 method for class 'shrinkr_prior_pred'
print(x, ...)
```

## Arguments

- x:

  A `shrinkr_prior_pred` object from
  [`sample_prior_predictive`](sample_prior_predictive.md).

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the input object `x`.

## See also

[`sample_prior_predictive`](sample_prior_predictive.md) for generating
samples, [`summary.shrinkr_prior_pred`](summary.shrinkr_prior_pred.md)
for detailed summaries

## Examples

``` r
priors <- list(
  mu = distributional::dist_normal(0, 5),
  tau = distributional::dist_truncated(distributional::dist_normal(0, 1), lower = 0)
)
prior_pred <- sample_prior_predictive(priors, n_groups = 3, n_draws = 50)
print(prior_pred)
#> == Prior Predictive Samples =========================
#> 
#> Draws:   50 
#> Groups:  3 
#> 
#> Prior specifications:
#>   mu:   N(0, 25) 
#>   tau:  N(0, 1)[0,Inf) 
#> 
#> Hyperparameter summaries:
#>   mu:  mean = -0.157 , sd = 5.474 , range = [ -14.445 , 12.488 ]
#>   tau: mean = 0.879 , sd = 0.658 , range = [ 0.01 , 2.875 ]
#> 
#> Theta (by group):
#>    group1 : mean = -0.209 , sd = 5.573 
#>    group2 : mean = -0.06 , sd = 5.638 
#>    group3 : mean = -0.019 , sd = 5.7 
#> 
#> -----------------------------------------------------
#> Use plot() to visualize
#> Use as.data.frame() for tidy format
#> Use summary() for detailed statistics
```
