# Summarize mu and tau hyperparameters

Computes posterior summaries for the hierarchical hyperparameters (`mu`,
`tau`, and `tau_squared`). Returns a data frame with one row per
parameter containing posterior means, standard deviations, quantiles,
and convergence diagnostics.

This is a focused alternative to `summary(fit)`, which returns summaries
for all parameters including theta.

## Usage

``` r
summarise_mu_tau(fit, probs = c(0.025, 0.5, 0.975), measures = NULL)

summarize_mu_tau(fit, probs = c(0.025, 0.5, 0.975), measures = NULL)
```

## Arguments

- fit:

  A `shrinkr_fit` object from [`shrink()`](shrink.md).

- probs:

  Numeric vector of quantiles to compute. Default is
  `c(0.025, 0.5, 0.975)` for 95% credible intervals.

- measures:

  Optional character vector or list of summary measures to compute. If
  `NULL`, uses mean, sd, and convergence diagnostics.

## Value

A data frame (tibble if available) with one row per parameter and
columns:

- `parameter`:

  Parameter name (`mu`, `tau`, or `tau_squared`)

- `mean`:

  Posterior mean

- `sd`:

  Posterior standard deviation

- `q2.5`, `q50`, `q97.5`:

  Quantiles (or custom quantiles from `probs`)

- `rhat`:

  R-hat convergence diagnostic

- `ess_bulk`:

  Effective sample size (bulk)

- `ess_tail`:

  Effective sample size (tail)

## See also

[`shrink()`](shrink.md) for fitting models,
[`extract_mu_tau()`](extract_mu_tau.md) for raw hyperparameter draws,
[`summarise_theta()`](summarise_theta.md) for group-level summaries

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- shrink(mixture = mix, hierarchical_priors = priors)

# Basic summary
mu_tau_summary <- summarise_mu_tau(fit)
print(mu_tau_summary)

# Custom quantiles
mu_tau_summary <- summarise_mu_tau(fit, probs = c(0.05, 0.5, 0.95))

# Custom measures
mu_tau_summary <- summarise_mu_tau(fit, measures = c("mean", "median", "mad"))
} # }
```
