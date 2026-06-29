# Summarize theta parameters by group

Computes posterior summaries for subgroup effects (theta parameters).
Returns a data frame with one row per group containing posterior means,
standard deviations, quantiles, and convergence diagnostics.

This is a focused alternative to `summary(fit)`, which returns summaries
for all parameters including mu and tau.

## Usage

``` r
summarise_theta(
  fit,
  probs = c(0.025, 0.5, 0.975),
  group_names = NULL,
  measures = NULL
)

summarize_theta(
  fit,
  probs = c(0.025, 0.5, 0.975),
  group_names = NULL,
  measures = NULL
)
```

## Arguments

- fit:

  A `shrinkr_fit` object from [`shrink()`](shrink.md).

- probs:

  Numeric vector of quantiles to compute. Default is
  `c(0.025, 0.5, 0.975)` for 95% credible intervals.

- group_names:

  Optional character vector of length G to label groups. If `NULL`, uses
  names from `fit$data$vars` or defaults to "group1", etc.

- measures:

  Optional character vector or list of summary measures to compute. If
  `NULL`, uses mean, sd, and convergence diagnostics.

## Value

A data frame (tibble if available) with one row per group and columns:

- `group`:

  Group identifier

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
[`summarise_mu_tau()`](summarise_mu_tau.md) for hyperparameter
summaries, [`theta_contrasts()`](theta_contrasts.md) for computing
contrasts

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
summarise_theta(fit)
#> # A tibble: 3 × 6
#>   group    mean     sd   q2.5    q50 q97.5
#>   <chr>   <dbl>  <dbl>  <dbl>  <dbl> <dbl>
#> 1 group1 0.0139 0.0810 -0.109 0.0114 0.172
#> 2 group2 0.310  0.105   0.146 0.305  0.529
#> 3 group3 0.512  0.0911  0.359 0.530  0.641
summarise_theta(fit, probs = c(0.05, 0.5, 0.95))
#> # A tibble: 3 × 6
#>   group    mean     sd     q5    q50   q95
#>   <chr>   <dbl>  <dbl>  <dbl>  <dbl> <dbl>
#> 1 group1 0.0139 0.0810 -0.105 0.0114 0.146
#> 2 group2 0.310  0.105   0.172 0.305  0.518
#> 3 group3 0.512  0.0911  0.371 0.530  0.623
summarise_theta(fit, group_names = c("Control", "A", "B"))
#> # A tibble: 3 × 6
#>   group     mean     sd   q2.5    q50 q97.5
#>   <chr>    <dbl>  <dbl>  <dbl>  <dbl> <dbl>
#> 1 Control 0.0139 0.0810 -0.109 0.0114 0.172
#> 2 A       0.310  0.105   0.146 0.305  0.529
#> 3 B       0.512  0.0911  0.359 0.530  0.641
```
