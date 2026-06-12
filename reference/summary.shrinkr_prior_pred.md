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
if (FALSE) { # \dontrun{
library(distributional)

priors <- list(
  mu = dist_normal(0, 5),
  tau = dist_truncated(dist_normal(0, 1), lower = 0)
)

prior_pred <- sample_prior_predictive(priors, n_groups = 4, n_draws = 1000)

# Get detailed summary
summ <- summary(prior_pred)

# View hyperparameter summaries
summ$hyperparameters

# View theta summaries
summ$theta

# Custom quantiles
summary(prior_pred, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
} # }
```
