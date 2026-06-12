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
if (FALSE) { # \dontrun{
fit <- shrink(mixture = mix, hierarchical_priors = priors)

# Pairwise contrast: group 2 vs group 1
L <- matrix(c(-1, 1, 0, 0), nrow = 1)
contrast <- theta_contrasts(fit, L, labels = "Trt2_vs_Trt1")

# Multiple contrasts
L <- rbind(
  c(-1, 1, 0, 0),   # Group 2 vs 1
  c(-1, 0, 1, 0),   # Group 3 vs 1
  c(0, -1, 1, 0)    # Group 3 vs 2
)
contrasts <- theta_contrasts(
  fit, L,
  labels = c("2vs1", "3vs1", "3vs2")
)

# Summarize contrasts
library(posterior)
summarise_draws(contrasts)

# Visualize
library(bayesplot)
mcmc_areas(contrasts, regex_pars = ".*")

# Probability that contrast > 0
mean(contrasts$`2vs1` > 0)
} # }
```
