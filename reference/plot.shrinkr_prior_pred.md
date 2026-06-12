# Plot prior predictive samples

Visualizes the prior predictive distribution for hyperparameters (mu and
tau) and subgroup effects (theta). Requires the ggplot2 package.

## Usage

``` r
# S3 method for class 'shrinkr_prior_pred'
plot(x, type = c("both", "hyperparameters", "theta"), ...)
```

## Arguments

- x:

  A `shrinkr_prior_pred` object from
  [`sample_prior_predictive`](sample_prior_predictive.md).

- type:

  Character; type of plot. Options:

  - `"both"` - Both hyperparameter and theta plots (default)

  - `"hyperparameters"` - Density plots for mu and tau only

  - `"theta"` - Violin plots for theta by group only

- ...:

  Additional arguments (currently unused).

## Value

A ggplot2 object, or a list of two ggplot2 objects if `type = "both"`
and patchwork package is not available. If patchwork is available and
`type = "both"`, returns a combined plot.

## See also

[`sample_prior_predictive`](sample_prior_predictive.md) for generating
samples,
[`as.data.frame.shrinkr_prior_pred`](as.data.frame.shrinkr_prior_pred.md)
for extracting data

## Examples

``` r
if (FALSE) { # \dontrun{
library(distributional)

priors <- list(
  mu = dist_normal(0, 5),
  tau = dist_truncated(dist_student_t(3, 0, 1), lower = 0)
)

prior_pred <- sample_prior_predictive(priors, n_groups = 4, n_draws = 1000)

# Plot both hyperparameters and theta
plot(prior_pred)

# Plot only hyperparameters
plot(prior_pred, type = "hyperparameters")

# Plot only theta
plot(prior_pred, type = "theta")
} # }
```
