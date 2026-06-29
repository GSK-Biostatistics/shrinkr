# Linear combinations of theta

Computes posterior draws for linear combinations of subgroup effects.
Useful for pairwise contrasts (e.g., treatment vs control), weighted
averages, or any custom linear estimand involving theta parameters.

## Usage

``` r
theta_contrasts(fit, contrast_matrix, labels = NULL)
```

## Arguments

- fit:

  A `shrinkr_fit` object from [`shrink()`](shrink.md).

- contrast_matrix:

  A numeric matrix L with `ncol(L) = G` (number of groups) and
  `nrow(L) = M` (number of contrasts). Each row defines one linear
  combination: \$\$contrast_i = L\_{i1}\theta_1 + L\_{i2}\theta_2 +
  \ldots + L\_{iG}\theta_G\$\$

- labels:

  Optional character vector of length M to name the contrasts. If
  `NULL`, uses "contrast1", "contrast2", etc.

## Value

A
[`posterior::draws_df`](https://mc-stan.org/posterior/reference/draws_df.html)
with columns `.chain`, `.iteration`, `.draw`, and one column per
contrast.

## See also

[`shrink()`](shrink.md) for fitting models,
[`summarise_theta()`](summarise_theta.md) for basic theta summaries

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
L <- matrix(c(-1, 1, 0), nrow = 1)
contrast <- theta_contrasts(fit, L, labels = "group2_vs_group1")
posterior::summarise_draws(contrast)
#> # A tibble: 1 × 10
#>   variable         mean median    sd    mad     q5   q95  rhat ess_bulk ess_tail
#>   <chr>           <dbl>  <dbl> <dbl>  <dbl>  <dbl> <dbl> <dbl>    <dbl>    <dbl>
#> 1 group2_vs_grou… 0.296  0.304 0.118 0.0756 0.0805 0.436  1.06     17.1     20.4
```
