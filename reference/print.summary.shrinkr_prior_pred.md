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
if (FALSE) { # \dontrun{
library(distributional)

priors <- list(
  mu = dist_normal(0, 5),
  tau = dist_truncated(dist_normal(0, 1), lower = 0)
)

prior_pred <- sample_prior_predictive(priors, n_groups = 4, n_draws = 1000)
summ <- summary(prior_pred)

# Print summary
print(summ)

# With more decimal places
print(summ, digits = 4)
} # }
```
