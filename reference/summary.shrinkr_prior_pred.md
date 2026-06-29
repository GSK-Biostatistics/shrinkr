# Summary statistics for prior predictive samples

Computes comprehensive summary statistics for hyperparameters (mu and
tau) and theta parameters from prior predictive samples.

## Usage

``` r
# S3 method for class 'shrinkr_prior_pred'
summary(object, probs = c(0.025, 0.5, 0.975), ...)
```

## Arguments

- object:

  A `shrinkr_prior_pred` object from
  [`sample_prior_predictive`](sample_prior_predictive.md).

- probs:

  Numeric vector of quantiles to compute. Default is c(0.025, 0.5,
  0.975) for 95% credible intervals.

- ...:

  Additional arguments (currently unused).

## Value

A list with class "summary.shrinkr_prior_pred" containing:

- hyperparameters:

  Data frame with summary statistics for mu and tau

- theta:

  Data frame with summary statistics for each group's theta

Each data frame includes columns for mean, sd, and the requested
quantiles.

## See also

[`sample_prior_predictive`](sample_prior_predictive.md) for generating
samples, [`print.shrinkr_prior_pred`](print.shrinkr_prior_pred.md) for
quick overview

## Examples

``` r
priors <- list(
  mu = distributional::dist_normal(0, 5),
  tau = distributional::dist_truncated(distributional::dist_normal(0, 1), lower = 0)
)
prior_pred <- sample_prior_predictive(priors, n_groups = 3, n_draws = 50)
summ <- summary(prior_pred)
summ$theta
#> # A tibble: 3 × 6
#>   group     mean    sd  q2.5   q50.0 q97.5
#>   <chr>    <dbl> <dbl> <dbl>   <dbl> <dbl>
#> 1 group1 -0.209   5.57 -10.7 -0.341  10.8 
#> 2 group2 -0.0604  5.64 -10.9  0.190  10.0 
#> 3 group3 -0.0191  5.70 -10.7 -0.0812  9.51
```
