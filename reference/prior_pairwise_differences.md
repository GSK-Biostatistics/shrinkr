# Compute prior predictive pairwise differences \|theta_i - theta_j\|

Computes the prior-implied distribution of absolute pairwise differences
between subgroup effects. This is useful for calibrating priors: if your
prior implies that subgroup differences of 5 units are common but
clinical relevance starts at 0.5, your prior may be too diffuse.

This implements the recommendation from the SDSIH Vignettes Library:
inspect the prior distribution of \\\|\theta_i - \theta_j\|\\ when
choosing priors for Bayesian hierarchical models.

## Usage

``` r
prior_pairwise_differences(prior_pred)
```

## Arguments

- prior_pred:

  A `shrinkr_prior_pred` object from
  [`sample_prior_predictive`](sample_prior_predictive.md).

## Value

A list with class `"shrinkr_prior_contrasts"` containing:

- differences:

  Data frame with columns `pair`, `abs_diff`, and `.draw`

- summary:

  Data frame with per-pair summary statistics

- overall_summary:

  Named numeric vector of quantiles across all pairs

- n_pairs:

  Number of unique pairs

- n_draws:

  Number of prior predictive draws

- group_names:

  Group labels used

## See also

[`sample_prior_predictive`](sample_prior_predictive.md) for generating
prior predictive samples,
[`plot.shrinkr_prior_contrasts`](plot.shrinkr_prior_contrasts.md) for
visualizing the result

## Examples

``` r
priors <- list(
  mu = distributional::dist_normal(0, 5),
  tau = distributional::dist_truncated(distributional::dist_normal(0, 1), lower = 0)
)
prior_pred <- sample_prior_predictive(priors, n_groups = 3, n_draws = 50)
pw <- prior_pairwise_differences(prior_pred)
print(pw)
#> == Prior Predictive: Pairwise |theta_i - theta_j| ==
#> 
#> Groups:  3 
#> Pairs:   3 
#> Draws:   50 
#> 
#> Overall quantiles of |theta_i - theta_j|:
#>    q2.5 = 0.007, q25 = 0.159, q50 = 0.448, q75 = 1.22, q97.5 = 5.372 
#> 
#> Per-pair summary:
#> # A tibble: 3 × 6
#>   pair             median    q2.5 q97.5 prob_gt_0.5 prob_gt_1
#>   <chr>             <dbl>   <dbl> <dbl>       <dbl>     <dbl>
#> 1 group1 vs group2  0.391 0.00796  5.51        0.48      0.36
#> 2 group1 vs group3  0.474 0.0139   8.97        0.48      0.32
#> 3 group2 vs group3  0.463 0.00648  3.77        0.46      0.3 
#> 
#> -----------------------------------------------------
#> Use plot() to visualize
pw$overall_summary
#>        q2.5         q25         q50         q75       q97.5 
#> 0.006571204 0.158913097 0.447961911 1.220236994 5.372119638 
```
