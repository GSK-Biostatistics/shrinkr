# Convert prior predictive samples to data frame

Converts prior predictive samples to a tidy long-format data frame
suitable for analysis and visualization with tidyverse tools.

## Usage

``` r
# S3 method for class 'shrinkr_prior_pred'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `shrinkr_prior_pred` object from
  [`sample_prior_predictive`](sample_prior_predictive.md).

- row.names:

  Ignored (for S3 consistency).

- optional:

  Ignored (for S3 consistency).

- ...:

  Additional arguments (currently unused).

## Value

A data frame (or tibble if tibble package is available) with columns:

- .draw:

  Draw number (1 to n_draws)

- group:

  Group name

- theta:

  Sampled group-level effect

- mu:

  Sampled global mean for this draw

- tau:

  Sampled heterogeneity parameter for this draw

## See also

[`sample_prior_predictive`](sample_prior_predictive.md) for generating
prior predictive samples

## Examples

``` r
if (FALSE) { # \dontrun{
library(distributional)

priors <- list(
  mu = dist_normal(0, 5),
  tau = dist_truncated(dist_normal(0, 1), lower = 0)
)

prior_pred <- sample_prior_predictive(priors, n_groups = 3, n_draws = 500)

# Convert to data frame
df <- as.data.frame(prior_pred)
head(df)

# Use with dplyr
library(dplyr)
df %>%
  group_by(group) %>%
  summarise(mean_theta = mean(theta), sd_theta = sd(theta))
} # }
```
