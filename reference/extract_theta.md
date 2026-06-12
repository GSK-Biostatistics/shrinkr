# Extract theta (group-level effect) parameters

Extracts posterior draws for the group-level effects (theta parameters)
from a fitted shrinkage model. This is the hierarchically shrunk version
of the subgroup effects.

## Usage

``` r
extract_theta(x, ...)
```

## Arguments

- x:

  A `shrinkr_fit` object from [`shrink()`](shrink.md).

- ...:

  Additional arguments passed to
  [`as_draws_df.shrinkr_fit()`](as_draws_df.shrinkr_fit.md).

## Value

A
[`posterior::draws_df`](https://mc-stan.org/posterior/reference/draws_df.html)
with columns:

- `.chain`:

  Chain index

- `.iteration`:

  Iteration within chain

- `.draw`:

  Overall draw index

- `theta[1]`, `theta[2]`, ...:

  Group-level effects

## See also

[`shrink()`](shrink.md) for fitting models,
[`extract_mu_tau()`](extract_mu_tau.md) for hyperparameters,
[`summarise_theta()`](summarise_theta.md) for summary statistics,
[`theta_contrasts()`](theta_contrasts.md) for pairwise comparisons

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- shrink(mixture = mix, hierarchical_priors = priors)

# Extract theta draws
theta_draws <- extract_theta(fit)

# Summarize
library(posterior)
summarise_draws(theta_draws)

# Visualize
library(bayesplot)
mcmc_intervals(theta_draws)

# Compare to summaries
theta_summary <- summarise_theta(fit)
print(theta_summary)
} # }
```
