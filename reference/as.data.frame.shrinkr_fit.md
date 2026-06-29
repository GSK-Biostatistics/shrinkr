# Convert shrinkr_fit to data.frame

Extracts posterior draws as a regular data frame. This is a convenience
wrapper around `as_draws_df()` that returns a plain data.frame.

## Usage

``` r
# S3 method for class 'shrinkr_fit'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  variables = NULL,
  include_internals = FALSE,
  ...
)
```

## Arguments

- x:

  A `shrinkr_fit` object from [`shrink()`](shrink.md).

- row.names:

  NULL or character vector giving row names.

- optional:

  Logical; if `TRUE`, setting row names and converting column names is
  optional.

- variables:

  Character vector of parameter names to extract. If `NULL`, returns all
  user-facing parameters (excludes internals).

- include_internals:

  Logical; if `TRUE`, includes internal Stan parameters. Default
  `FALSE`.

- ...:

  Additional arguments passed to `as_draws_df()`.

## Value

A data.frame with columns for chain, iteration, draw, and requested
parameters.

## See also

[`as_draws_df.shrinkr_fit()`](as_draws_df.shrinkr_fit.md) for posterior
package format, [`extract_mu_tau()`](extract_mu_tau.md) for
hyperparameters only

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
draws_df <- as.data.frame(fit)
head(draws_df)
#>          mu       tau    theta[1]  theta[2]  theta[3] tau_squared .chain
#> 1 0.1686773 0.3275693 -0.01645236 0.5401618 0.4431331  0.10730166      1
#> 2 0.2091822 0.3234641 -0.02533617 0.2960760 0.4864821  0.10462902      1
#> 3 0.1582186 0.3022369  0.06969634 0.3689739 0.6178087  0.09134717      1
#> 4 0.2797640 0.2403194  0.05566632 0.3028002 0.3476433  0.05775344      1
#> 5 0.2164754 0.3185948 -0.06887557 0.2256727 0.5593946  0.10150263      1
#> 6 0.1589766 0.2983161 -0.07074952 0.3188792 0.5332950  0.08899252      1
#>   .iteration .draw
#> 1          1     1
#> 2          2     2
#> 3          3     3
#> 4          4     4
#> 5          5     5
#> 6          6     6
mu_tau_df <- as.data.frame(fit, variables = c("mu", "tau"))
```
