# Convert prior predictive samples to data frame

Converts prior predictive samples to a tidy long-format data frame
suitable for analysis and visualization with tidyverse tools.

## Usage

``` r
# S3 method for class 'shrinkr_prior_pred'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `shrinkr_prior_pred` object from
  [`sample_prior_predictive`](sample_prior_predictive.md).

- row.names:

  Ignored (for S3 consistency).

- optional:

  Ignored (for S3 consistency).

- ...:

  Additional arguments (currently unused).

## Value

A data frame (or tibble if tibble package is available) with columns:

- .draw:

  Draw number (1 to n_draws)

- group:

  Group name

- theta:

  Sampled group-level effect

- mu:

  Sampled global mean for this draw

- tau:

  Sampled heterogeneity parameter for this draw

## See also

[`sample_prior_predictive`](sample_prior_predictive.md) for generating
prior predictive samples

## Examples

``` r
priors <- list(
  mu = distributional::dist_normal(0, 5),
  tau = distributional::dist_truncated(distributional::dist_normal(0, 1), lower = 0)
)
prior_pred <- sample_prior_predictive(priors, n_groups = 3, n_draws = 50)
df <- as.data.frame(prior_pred)
head(df)
#> # A tibble: 6 × 5
#>   .draw group  theta    mu   tau
#>   <int> <chr>  <dbl> <dbl> <dbl>
#> 1     1 group1 2.91   2.05 0.762
#> 2     1 group2 0.302  2.05 0.762
#> 3     1 group3 2.61   2.05 0.762
#> 4     2 group1 7.11   8.44 1.01 
#> 5     2 group2 9.37   8.44 1.01 
#> 6     2 group3 8.85   8.44 1.01 
```
