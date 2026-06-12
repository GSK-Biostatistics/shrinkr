# Plot prior predictive pairwise differences

Creates a density plot of \\\|\theta_i - \theta_j\|\\ from prior
predictive samples. Useful for calibrating hierarchical priors. Styling
matches [`plot.shrinkr_prior_pred()`](plot.shrinkr_prior_pred.md) for
visual consistency.

## Usage

``` r
# S3 method for class 'shrinkr_prior_contrasts'
plot(x, by_pair = FALSE, ...)
```

## Arguments

- x:

  A `shrinkr_prior_contrasts` object from
  [`prior_pairwise_differences`](prior_pairwise_differences.md).

- by_pair:

  Logical; if `TRUE`, facet by pair. If `FALSE` (default), pool all
  pairwise differences into a single plot.

- ...:

  Additional arguments (currently unused).

## Value

A ggplot2 object.

## See also

[`prior_pairwise_differences`](prior_pairwise_differences.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(distributional)

priors <- list(
  mu = dist_normal(0, 5),
  tau = dist_truncated(dist_normal(0, 1), lower = 0)
)

prior_pred <- sample_prior_predictive(priors, n_groups = 4, n_draws = 2000)
pw <- prior_pairwise_differences(prior_pred)

# Pooled density
plot(pw)

# Faceted by pair
plot(pw, by_pair = TRUE)
} # }
```
