# Summary method for shrinkr_fit

Comprehensive posterior summary including hyperparameters and all
subgroup effects with convergence diagnostics.

## Usage

``` r
# S3 method for class 'shrinkr_fit'
summary(
  object,
  probs = c(0.025, 0.5, 0.975),
  group_names = NULL,
  digits = 3,
  ...
)
```

## Arguments

- object:

  A `shrinkr_fit` object.

- probs:

  Numeric vector of quantiles to compute. Default is
  `c(0.025, 0.5, 0.975)`.

- group_names:

  Optional character vector to label groups.

- digits:

  Number of digits to display. Default is 3.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns a list with `mu_tau` and `theta` summary tables.
Prints formatted output.

## See also

[`shrink()`](shrink.md) for fitting models,
[`summarise_theta()`](summarise_theta.md) for theta-only summaries
