# Create a mixture prior

Creates a mixture distribution using
[`distributional::dist_mixture()`](https://pkg.mitchelloharawild.com/distributional/reference/dist_mixture.html).
All standard distributional operations (sampling, density, quantiles,
formatting) work automatically.

## Usage

``` r
prior_mixture(..., weights = NULL)
```

## Arguments

- ...:

  Component distributions (from distributional package)

- weights:

  Mixture weights (normalized automatically if not summing to 1)

## Value

A `distributional` mixture distribution object

## Examples

``` r
mix <- prior_mixture(
  distributional::dist_normal(0, 0.1),
  distributional::dist_normal(0, 1),
  weights = c(0.7, 0.3)
)
```
