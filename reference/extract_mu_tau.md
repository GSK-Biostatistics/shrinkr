# Extract mu and tau parameters

Extracts posterior draws for the hyperparameters mu (global mean) and
tau (heterogeneity standard deviation) from a fitted shrinkage model.

## Usage

``` r
extract_mu_tau(x, ...)
```

## Arguments

- x:

  A `shrinkr_fit` object from [`shrink()`](shrink.md).

- ...:

  Additional arguments (currently unused).

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

- `mu`:

  Global mean parameter

- `tau`:

  Heterogeneity parameter

- `tau_squared`:

  Variance (tau^2)

## See also

[`shrink()`](shrink.md) for fitting models,
[`summarise_mu_tau()`](summarise_mu_tau.md) for summary statistics,
[`as_draws_df.shrinkr_fit()`](as_draws_df.shrinkr_fit.md) for all
parameters

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- shrink(mixture = mix, hierarchical_priors = priors)

# Extract hyperparameter draws
mu_tau <- extract_mu_tau(fit)

# Summarize
summarise_mu_tau(fit)

# Visualize
library(bayesplot)
mcmc_pairs(mu_tau, pars = c("mu", "tau"))
} # }
```
