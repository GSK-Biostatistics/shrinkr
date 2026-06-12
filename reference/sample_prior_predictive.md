# Sample from prior predictive distribution

Generates samples from the prior predictive distribution for the
hierarchical shrinkage model. Useful for prior elicitation and
sensitivity analysis.

The generative process is:

1.  Sample mu from p(mu)

2.  Sample tau from p(tau)

3.  Sample theta_i ~ N(mu, tau) for each group i

## Usage

``` r
sample_prior_predictive(
  hierarchical_priors,
  n_groups,
  n_draws = 1000,
  group_names = NULL
)
```

## Arguments

- hierarchical_priors:

  Named list with `mu` and `tau` distributional objects from the
  `distributional` package.

- n_groups:

  Integer; number of subgroups (G).

- n_draws:

  Integer; number of prior predictive samples to draw. Default 1000.

- group_names:

  Optional character vector of length `n_groups` to label groups.

## Value

A list with class "shrinkr_prior_pred" containing:

- mu:

  Vector of mu draws

- tau:

  Vector of tau draws

- theta:

  Matrix of theta draws (n_draws x n_groups)

- implied_range:

  Vector of ranges (max - min) of theta across groups for each draw

- implied_sd:

  Vector of standard deviations of theta across groups for each draw

- group_names:

  Group labels

- n_draws:

  Number of draws

- n_groups:

  Number of groups

- priors:

  The hierarchical_priors specification used

## See also

[`shrink`](shrink.md) for fitting the hierarchical model,
[`plot.shrinkr_prior_pred`](plot.shrinkr_prior_pred.md) for visualizing
prior predictive samples

## Examples

``` r
if (FALSE) { # \dontrun{
library(distributional)

# Specify priors
priors <- list(
  mu = dist_normal(0, 5),
  tau = dist_truncated(dist_normal(0, 2.5), lower = 0)
)

# Sample from prior predictive
prior_pred <- sample_prior_predictive(
  hierarchical_priors = priors,
  n_groups = 4,
  n_draws = 1000
)

# Visualize
plot(prior_pred)

# Check implied spread of effects
cat("Median implied range:", median(prior_pred$implied_range), "\n")
cat("Median implied SD:", median(prior_pred$implied_sd), "\n")

# Extract as tidy data frame
prior_df <- as.data.frame(prior_pred)

# Check prior implied theta distribution
library(ggplot2)
ggplot(prior_df, aes(x = theta, y = group)) +
  geom_violin() +
  labs(title = "Prior predictive distribution for theta")
} # }
```
