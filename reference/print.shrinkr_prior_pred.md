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
if (FALSE) { # \dontrun{
library(distributional)

priors <- list(
  mu = dist_normal(0, 5),
  tau = dist_truncated(dist_normal(0, 1), lower = 0)
)

prior_pred <- sample_prior_predictive(priors, n_groups = 4, n_draws = 1000)

# Print summary
print(prior_pred)

# Or just type the object name
prior_pred
} # }
```
