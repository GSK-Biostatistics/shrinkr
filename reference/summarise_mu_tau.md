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
summarise_mu_tau(fit)
#> # A tibble: 3 × 6
#>   parameter     mean     sd   q2.5    q50 q97.5
#>   <chr>        <dbl>  <dbl>  <dbl>  <dbl> <dbl>
#> 1 mu          0.210  0.0457 0.122  0.218  0.278
#> 2 tau         0.300  0.0261 0.248  0.298  0.337
#> 3 tau_squared 0.0905 0.0153 0.0614 0.0890 0.114
summarise_mu_tau(fit, probs = c(0.05, 0.5, 0.95))
#> # A tibble: 3 × 6
#>   parameter     mean     sd     q5    q50   q95
#>   <chr>        <dbl>  <dbl>  <dbl>  <dbl> <dbl>
#> 1 mu          0.210  0.0457 0.155  0.218  0.276
#> 2 tau         0.300  0.0261 0.255  0.298  0.333
#> 3 tau_squared 0.0905 0.0153 0.0651 0.0890 0.111
```
