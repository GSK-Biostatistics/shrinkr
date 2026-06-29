# Plot shrinkage fit

Visualizes the hierarchical shrinkage model fit. Creates either:

- `"shrinkage"` — Shows pre/post-shrunk estimates with arrows

- `"diagnostics"` — Multi-panel view with hyperparameters and shrinkage
  factor

The shrinkage plot displays:

- **Pre-shrunk estimates** (hollow circles) — from stage 1 mixture or
  MLEs

- **Post-shrunk estimates** (filled circles) — posterior means of theta

- **Global mean** (dashed line) — posterior mean of mu

- **Credible intervals** (optional) — for shrunken estimates

- **Arrows** (optional) — showing direction and magnitude of shrinkage

## Usage

``` r
# S3 method for class 'shrinkr_fit'
plot(
  x,
  type = c("shrinkage", "diagnostics"),
  group_names = NULL,
  show_arrows = FALSE,
  show_intervals = TRUE,
  interval_prob = 0.95,
  point_size = 3,
  arrow_alpha = 0.6,
  dodge_width = 0.3,
  title = NULL,
  subtitle = NULL,
  ...
)
```

## Arguments

- x:

  A `shrinkr_fit` object from [`shrink()`](shrink.md).

- type:

  Character; type of plot. Options:

  - `"shrinkage"` - Basic shrinkage visualization (default)

  - `"diagnostics"` - Multi-panel with hyperparameters and shrinkage
    factor

- group_names:

  Optional character vector of length G to label groups. If `NULL`, uses
  names from `x$data$vars` or defaults to "group1", etc.

- show_arrows:

  Logical; draw arrows from pre-shrunk to post-shrunk estimates? Default
  `FALSE`. Only applies when `type = "shrinkage"`.

- show_intervals:

  Logical; show credible intervals for both pre-shrunk and post-shrunk
  estimates? Default `TRUE`. Only applies when `type = "shrinkage"`.

- interval_prob:

  Numeric; probability mass for credible intervals. Default 0.95 for 95%
  intervals. Only applies when `type = "shrinkage"`.

- point_size:

  Numeric; size of points. Default 3.

- arrow_alpha:

  Numeric; transparency of arrows (0-1). Default 0.6. Only applies when
  `show_arrows = TRUE`.

- dodge_width:

  Numeric; horizontal spacing between pre-shrunk and post-shrunk
  estimates in the side-by-side display. Default 0.3. Larger values
  increase separation between estimate types.

- title:

  Character; plot title. If `NULL`, uses default title.

- subtitle:

  Character; plot subtitle. If `NULL`, auto-generates from global mean
  and tau.

- ...:

  Additional arguments (currently unused).

## Value

A ggplot2 object (for `type = "shrinkage"`), or a patchwork object/list
(for `type = "diagnostics"`).

## See also

[`shrink()`](shrink.md) for fitting models,
[`extract_mu_tau()`](extract_mu_tau.md) for hyperparameter draws

## Examples

``` r
# Plotting requires a fitted shrinkr_fit object from shrink().
# The full example is not run because it fits a Stan model.
if (FALSE) { # \dontrun{
fit <- shrink(mixture = mix, hierarchical_priors = priors)
plot(fit)
plot(fit, show_arrows = TRUE, interval_prob = 0.95)
} # }
```
