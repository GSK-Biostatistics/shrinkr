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
set.seed(1)
samples <- list(
  group1 = matrix(rnorm(100, 0.0, 0.5), ncol = 1),
  group2 = matrix(rnorm(100, 0.5, 0.5), ncol = 1)
)
mix <- fit_mixture(samples, K_max = 2, verbose = FALSE)
df <- as.data.frame(mix)
head(df)
#> # A tibble: 2 × 5
#>   component variable weight   mean    sd
#>       <int> <chr>     <dbl>  <dbl> <dbl>
#> 1         1 group1        1 0.0544 0.462
#> 2         1 group2        1 0.481  0.462
```
