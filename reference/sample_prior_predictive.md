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
priors <- list(
  mu = distributional::dist_normal(0, 5),
  tau = distributional::dist_truncated(distributional::dist_normal(0, 1), lower = 0)
)
prior_pred <- sample_prior_predictive(priors, n_groups = 3, n_draws = 50)
median(prior_pred$implied_range)
#> [1] 0.9419793
head(as.data.frame(prior_pred))
#> # A tibble: 6 × 5
#>   .draw group   theta     mu   tau
#>   <int> <chr>   <dbl>  <dbl> <dbl>
#> 1     1 group1  0.326  0.305 0.305
#> 2     1 group2 -0.159  0.305 0.305
#> 3     1 group3  0.568  0.305 0.305
#> 4     2 group1 -3.65  -4.12  0.930
#> 5     2 group2 -4.45  -4.12  0.930
#> 6     2 group3 -4.57  -4.12  0.930
```
