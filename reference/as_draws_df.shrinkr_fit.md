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
set.seed(1)
draws <- data.frame(
  mu = rnorm(20, 0.2, 0.05),
  tau = abs(rnorm(20, 0.3, 0.03)),
  `theta[1]` = rnorm(20, 0.0, 0.1),
  `theta[2]` = rnorm(20, 0.3, 0.1),
  `theta[3]` = rnorm(20, 0.5, 0.1),
  check.names = FALSE
)
draws$tau_squared <- draws$tau^2
fit <- list(
  fit = posterior::as_draws_df(draws),
  data = list(
    G = 3, K = 1, centered = FALSE,
    vars = c("group1", "group2", "group3"),
    quantiles = data.frame(
      q2.5 = c(-0.20, 0.10, 0.30),
      q50 = c(0.00, 0.30, 0.50),
      q97.5 = c(0.20, 0.50, 0.70)
    )
  ),
  summary = posterior::summarise_draws(
    posterior::as_draws_df(draws),
    "mean", "sd",
    ~posterior::quantile2(., probs = c(0.025, 0.5, 0.975))
  ),
  diagnostics = list(n_divergent = 0, max_treedepth = 0, n_leapfrog = 0)
)
class(fit) <- "shrinkr_fit"
all_draws <- posterior::as_draws_df(fit)
posterior::variables(all_draws)
#> [1] "mu"          "tau"         "theta[1]"    "theta[2]"    "theta[3]"   
#> [6] "tau_squared"
theta_draws <- posterior::as_draws_df(fit, variables = c("theta[1]", "theta[2]"))
```
