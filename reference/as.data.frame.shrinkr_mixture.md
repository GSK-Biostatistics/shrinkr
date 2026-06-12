# Convert mixture fit to data frame

Converts a fitted mixture model to a tidy long-format data frame. This
is essentially an accessor for the `components` element.

## Usage

``` r
# S3 method for class 'shrinkr_mixture'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `shrinkr_mixture` object from [`fit_mixture`](fit_mixture.md).

- row.names:

  Ignored (for S3 consistency).

- optional:

  Ignored (for S3 consistency).

- ...:

  Additional arguments (currently unused).

## Value

A data frame (or tibble if available) containing the component
specifications with columns:

- component:

  Component number (1 to K)

- variable:

  Variable name

- weight:

  Component weight (mixing proportion)

- mean:

  Component mean for this variable

- sd:

  Component marginal standard deviation for this variable

## See also

[`fit_mixture`](fit_mixture.md) for fitting mixture models

## Examples

``` r
if (FALSE) { # \dontrun{
samples <- list(
  matrix(rnorm(2000, 0.0, 0.5), ncol = 1),
  matrix(rnorm(2000, 0.5, 0.5), ncol = 1),
  matrix(rnorm(2000, 1.0, 0.5), ncol = 1)
)

mix <- fit_mixture(samples, K_max = 5)

# Convert to data frame
df <- as.data.frame(mix)
head(df)

# Use with dplyr
library(dplyr)
df %>%
  group_by(variable) %>%
  summarise(weighted_mean = sum(weight * mean))
} # }
```
