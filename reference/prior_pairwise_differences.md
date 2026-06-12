# Compute prior predictive pairwise differences \|theta_i - theta_j\|

Computes the prior-implied distribution of absolute pairwise differences
between subgroup effects. This is useful for calibrating priors: if your
prior implies that subgroup differences of 5 units are common but
clinical relevance starts at 0.5, your prior may be too diffuse.

This implements the recommendation from the SDSIH Vignettes Library:
inspect the prior distribution of \\\|\theta_i - \theta_j\|\\ when
choosing priors for Bayesian hierarchical models.

## Usage

``` r
prior_pairwise_differences(prior_pred)
```

## Arguments

- prior_pred:

  A `shrinkr_prior_pred` object from
  [`sample_prior_predictive`](sample_prior_predictive.md).

## Value

A list with class `"shrinkr_prior_contrasts"` containing:

- differences:

  Data frame with columns `pair`, `abs_diff`, and `.draw`

- summary:

  Data frame with per-pair summary statistics

- overall_summary:

  Named numeric vector of quantiles across all pairs

- n_pairs:

  Number of unique pairs

- n_draws:

  Number of prior predictive draws

- group_names:

  Group labels used

## See also

[`sample_prior_predictive`](sample_prior_predictive.md) for generating
prior predictive samples,
[`plot.shrinkr_prior_contrasts`](plot.shrinkr_prior_contrasts.md) for
visualizing the result

## Examples

``` r
if (FALSE) { # \dontrun{
library(distributional)

priors <- list(
  mu = dist_normal(0, 5),
  tau = dist_truncated(dist_normal(0, 1), lower = 0)
)

prior_pred <- sample_prior_predictive(priors, n_groups = 4, n_draws = 2000)

# Compute pairwise differences
pw <- prior_pairwise_differences(prior_pred)
print(pw)
plot(pw)

# What range of pairwise differences does this prior imply?
pw$overall_summary
} # }
```
