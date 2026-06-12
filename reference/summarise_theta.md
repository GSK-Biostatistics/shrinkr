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
if (FALSE) { # \dontrun{
fit <- shrink(mixture = mix, hierarchical_priors = priors)

# Basic summary
theta_summary <- summarise_theta(fit)
print(theta_summary)

# Custom quantiles
theta_summary <- summarise_theta(fit, probs = c(0.05, 0.5, 0.95))

# With custom group names
theta_summary <- summarise_theta(
  fit, 
  group_names = c("Control", "Treatment A", "Treatment B")
)

# Custom measures
theta_summary <- summarise_theta(fit, measures = c("mean", "median", "mad"))
} # }
```
