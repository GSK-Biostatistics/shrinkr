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
theta_draws <- extract_theta(fit)
posterior::summarise_draws(theta_draws)
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> # A tibble: 6 × 10
#>   variable      mean median     sd    mad      q5   q95  rhat ess_bulk ess_tail
#>   <chr>        <dbl>  <dbl>  <dbl>  <dbl>   <dbl> <dbl> <dbl>    <dbl>    <dbl>
#> 1 mu          0.210  0.218  0.0457 0.0388  0.155  0.276 0.968     26.0     25.6
#> 2 tau         0.300  0.298  0.0261 0.0204  0.255  0.333 0.971     17.5     20.4
#> 3 theta[1]    0.0139 0.0114 0.0810 0.0789 -0.105  0.146 0.995     26.0     25.6
#> 4 theta[2]    0.310  0.305  0.105  0.0888  0.172  0.518 1.03      19.7     20.4
#> 5 theta[3]    0.512  0.530  0.0911 0.120   0.371  0.623 1.08      20.1     25.6
#> 6 tau_squared 0.0905 0.0890 0.0153 0.0124  0.0651 0.111 0.971     17.5     20.4
summarise_theta(fit)
#> # A tibble: 3 × 6
#>   group    mean     sd   q2.5    q50 q97.5
#>   <chr>   <dbl>  <dbl>  <dbl>  <dbl> <dbl>
#> 1 group1 0.0139 0.0810 -0.109 0.0114 0.172
#> 2 group2 0.310  0.105   0.146 0.305  0.529
#> 3 group3 0.512  0.0911  0.359 0.530  0.641
```
