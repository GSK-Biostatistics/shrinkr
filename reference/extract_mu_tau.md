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
mu_tau <- extract_mu_tau(fit)
summarise_mu_tau(fit)
#> # A tibble: 3 × 6
#>   parameter     mean     sd   q2.5    q50 q97.5
#>   <chr>        <dbl>  <dbl>  <dbl>  <dbl> <dbl>
#> 1 mu          0.210  0.0457 0.122  0.218  0.278
#> 2 tau         0.300  0.0261 0.248  0.298  0.337
#> 3 tau_squared 0.0905 0.0153 0.0614 0.0890 0.114
posterior::summarise_draws(mu_tau)
#> Warning: The ESS has been capped to avoid unstable estimates.
#> # A tibble: 3 × 10
#>   variable      mean median     sd    mad     q5   q95  rhat ess_bulk ess_tail
#>   <chr>        <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl>    <dbl>    <dbl>
#> 1 mu          0.210  0.218  0.0457 0.0388 0.155  0.276 0.968     26.0     25.6
#> 2 tau         0.300  0.298  0.0261 0.0204 0.255  0.333 0.971     17.5     20.4
#> 3 tau_squared 0.0905 0.0890 0.0153 0.0124 0.0651 0.111 0.971     17.5     20.4
```
