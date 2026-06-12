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
if (FALSE) { # \dontrun{
fit <- shrink(mixture = mix, hierarchical_priors = priors)

# Extract as data frame
draws_df <- as.data.frame(fit)
head(draws_df)

# Just mu and tau
mu_tau_df <- as.data.frame(fit, variables = c("mu", "tau"))
} # }
```
