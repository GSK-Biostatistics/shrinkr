# Convert shrinkr_fit to draws_df

Extracts posterior draws in tidy format using the `posterior` package.
By default returns user-facing parameters (mu, tau, theta, etc.) and
excludes internal parameterization details. Set
`include_internals = TRUE` to access all parameters including theta_c
and z.

## Usage

``` r
# S3 method for class 'shrinkr_fit'
as_draws_df(x, variables = NULL, include_internals = FALSE, ...)
```

## Arguments

- x:

  A `shrinkr_fit` object from [`shrink()`](shrink.md).

- variables:

  Character vector of parameter names to extract. Options include:

  - `"mu"` - Global mean

  - `"tau"` - Heterogeneity SD

  - `"tau_squared"` - Heterogeneity variance

  - `"theta"` or `"theta[i]"` - Subgroup effects

  If `NULL` (default), returns all user-facing parameters.

- include_internals:

  Logical; if `TRUE`, includes internal Stan parameters (`theta_c`, `z`)
  used for parameterization. Default `FALSE`. Only applies when
  `variables = NULL`.

- ...:

  Additional arguments passed to
  [`posterior::as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html).

## Value

A
[`posterior::draws_df`](https://mc-stan.org/posterior/reference/draws_df.html)
with columns for chain, iteration, draw, and requested parameters.

## See also

[`shrink()`](shrink.md) for fitting models,
[`extract_mu_tau()`](extract_mu_tau.md) for hyperparameters only

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- shrink(mixture = mix, hierarchical_priors = priors)

# User-facing parameters only (default)
all_draws <- as_draws_df(fit)
variables(all_draws)  # mu, tau, theta[1], ..., tau_squared

# Include internal parameters for diagnostics
all_draws_internal <- as_draws_df(fit, include_internals = TRUE)
variables(all_draws_internal)  # includes theta_c, z

# Just theta parameters
theta_draws <- as_draws_df(fit, variables = "theta")

# Specific thetas
theta12_draws <- as_draws_df(fit, variables = c("theta[1]", "theta[2]"))

# Work with draws
library(posterior)
summarise_draws(all_draws)
} # }
```
